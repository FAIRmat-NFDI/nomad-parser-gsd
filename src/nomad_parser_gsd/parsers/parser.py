import os
from typing import TYPE_CHECKING, Iterable, Union

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )
    from nomad.metainfo import (
        Context,
        Section,
    )
    from structlog.stdlib import (
        BoundLogger,
    )

from ase.symbols import symbols2numbers
from ase.utils import formula_hill
from collections import defaultdict
import nomad_simulations.schema_packages.properties.energies as energy_module
import nomad_simulations.schema_packages.properties.forces as force_module
import numpy as np

# from nomad.parsing.parser import MatchingParser
import structlog
from atomisticparsers.utils import MDParser
from nomad.config import config
from nomad.datamodel.metainfo.workflow import Workflow
from nomad.parsing.file_parser import FileParser
from nomad.units import ureg
from nomad_simulations.schema_packages.general import (
    Program,
    Simulation,
)

# nomad-simulations
from nomad_simulations.schema_packages.general import Program as BaseProgram
from nomad_simulations.schema_packages.outputs import TotalEnergy, TotalForce

# # nomad-parser-gsd
# from nomad_parser_gsd.schema_packages.schema import (
#     Author,
#     EnergyEntry,
#     ForceEntry,
#     ModelSystem,
#     OutputsEntry,
#     ParamEntry,
#     Stress,
#     TrajectoryOutputs,
# )

logging = structlog.get_logger()
try:
    import gsd.fl as gsdfl
    import gsd.hoomd as gsdhoomd
    from gsd.hoomd import HOOMDTrajectory
    # import gsd.pygsd as gsdpy
except ImportError:
    logging.warn('Required module gsd.hoomd not found.')
    gsdhoomd = False
    # gsdpy = False

configuration = config.get_plugin_entry_point(
    'nomad_parser_gsd.parsers:parser_entry_point'
)

MOL = 6.022140857e23


class GSDFileParser(FileParser):
    def __init__(self):
        super().__init__(None)

    @property
    def filegsd(self):
        if self._file_handler is None:
            try:
                self._file_handler = gsdhoomd.open(name=self.mainfile, mode='r')
            except Exception:
                self.logger.error('Error reading gsd file.')

            if type(self._file_handler) is not HOOMDTrajectory:
                self.logger.error(
                    'Uknown GSD file object, only HOOMDTrajectory objects are supported.'
                )
        return self._file_handler

    def get_value(self, group, path: str, default=None):
        """
        Extracts group or dataset from group object based on path, and returns default if not defined.
        """
        section_segments = path.split('.')
        for section in section_segments:
            try:
                value = getattr(group, section)
                group = value
            except AttributeError:
                return

        return value if value is not None else default

    def parse(self, path: str = None, **kwargs):
        frame_path = '.'.join(path.split('.')[1:])
        frame = kwargs.get('frame', None)

        if not frame:
            value = None
        else:
            value = self.get_value(frame, frame_path)

        self._results[path] = value


class GSDParser(MDParser):
    def __init__(self, *args, **kwargs):
        super().__init__()
        self._data_parser = GSDFileParser()
        self._maindir = None
        self._gsd_files = None
        self._basename = None
        self._program_dict = None
        self._n_frames = None
        self._first_frame = True
        self._time_unit = ureg.picosecond
        # ['N', 'position', 'orientation', 'types', 'typeid', 'mass', 'charge',
        # 'diameter', 'body', 'moment_inertia', 'velocity', 'angmom', 'image',
        # 'type_shapes']
        self._nomad_to_particles_group_map = {
            'n_atoms': 'N',
            'positions': 'position',
            'velocities': 'velocity',
            'forces': None,
            'labels': 'types',
            'label': None,
            'mass': 'mass',
            'charge': 'charge',
        }
        self._nomad_to_box_group_map = {
            'step': 'step',
            'lattice_vectors': '_box',
            'periodic_boundary_conditions': None,
            'dimensionality': 'dimensions',
        }

    # Load GSD file as file layer object to access generating program name and version.
    def get_program_info(self):
        _file_layer = gsdfl.open(name=self.mainfile, mode='r')
        _application = _file_layer.application.split(' ')

        _program_info = dict()

        _program_info['name'] = _application[0]
        _program_info['version'] = _application[1]
        _program_info['schema'] = (
            _file_layer.schema
        )  # Name of the data schema. In example: 'hoomd'
        _program_info['schema_version'] = (
            _file_layer.schema_version
        )  # TODO: tuple(major, minor). If keep, format to major.minor?
        _program_info['gsd_version'] = (
            _file_layer.gsd_version
        )  # TODO: tuple(major, minor). Format?

        self._n_frames = _file_layer.nframes  # ? Not ideal to access in get_program, but only file layer object has attribute nframes

        _file_layer.close()

        return _program_info

    """
    Logged data encompasses values computed at simulation time that are too expensive 
    or cumbersome to re-compute in post processing. This specification does not define 
    specific chunk names or define logged data. Users may select any valid name for 
    logged data chunks as appropriate for their workflow.
    """

    # Additional data are stored in the log dictionary as numpy arrays:
    def get_logged_info(self):
        try:
            return gsdhoomd.read_log(name=self.mainfile, scalar_only=False)
        except FileNotFoundError:
            self.logger.warning(
                'No additional logged data found, no user-defined data will be stored.'
            )
            return dict()

    # ? Will connectivity information still be used in this way?
    # connectivity
    # \-- particles_group  #! Entire topology
    #     \-- <group_1>
    #       |    \-- (type): String[]  #! molecule_group, molecule, monomer_group, monomer
    #       |    \-- (formula): String[]
    #       |    \-- indices: Integer[]  #! list of integer indices corresponding to all particles belonging to this group
    #       |    \-- (is_molecule): Bool
    #  |    |    \-- (<custom_dataset>): <type>[]
    #       |    \-- (particles_group): #! subgroups must be a subset of the grouping at the previous level of the hierarchy
    #       |        \-- ...  #! The particles_group hierarchy ends at the level of individual particles
    #! (i.e., individual particles are not stored, since this information is already contained within the particles group).
    #       \-- <group_2>
    #           \-- ...
    def get_particles_group(self, connectivity):
        # _interaction_keys = ['M', 'N', 'types', 'typeid', 'group']  # '_default_value'
        # _groups = defaultdict(list)
        # for idx, typeid in enumerate(self._particle_data_dict['typeid']):
        #     _groups[self._particle_data_dict['types'][typeid]].append(idx)

        _graph = defaultdict(list)
        for bond in connectivity['bonds']:
            _graph[bond[0]].append(bond[1])
            _graph[bond[1]].append(bond[0])

        # Depth-first search for connected components
        # -> extract molecules from bond list
        def dfs(node, visited, component):
            visited.add(node)
            component.append(node)
            for neighbor in _graph[node]:
                if neighbor not in visited:
                    dfs(neighbor, visited, component)

        _visited = set()
        _components = []
        # Traverse each node in the graph
        for node in _graph:
            if node not in _visited:
                _component = []
                dfs(node, _visited, _component)
                _components.append(_component)

        #! individual particles are not stored, since this information is already contained within the particles group
        # # Get all particles that are not part of any bond
        # _not_bond = [
        #     idx
        #     for idx, _ in enumerate(self._particle_data_dict['typeid'])
        #     if idx not in _visited
        # ]

        _particles_group = dict()
        _group_idx = 0

        def add_groups(_group_idx, components):  # _components, _not_bond
            _all_types = dict()
            _idx_dict = defaultdict(list)
            for component in components:
                _types = list()
                for idx in component:
                    _types.append(self._particle_data_dict['typeid'][idx])
                _types = tuple(sorted(_types))

                if _types not in set(_all_types.values()):
                    _all_types[_types] = f'group{_group_idx}'
                    _idx_dict[f'group{_group_idx}'].extend(component)
                    _group_idx += 1
                else:
                    group_name = _all_types[_types]
                    _idx_dict[group_name].extend(component)
            print(len(_components), len(_idx_dict.keys()))
            # _particles_group[f'group{_group_idx}'] = _group

        add_groups(_group_idx, _components)

        return _particles_group

    # connectivity
    # \-- (bonds): Integer[N_part][2] #! Need to be lists of tuples
    # \-- (angles): Integer[N_part][3]
    # \-- (dihedrals): Integer[N_part][4]
    # \-- (impropers): Integer[N_part][4]
    # \-- (<custom_interaction>): Integer[N_part][m]
    # \-- (particles_group)
    #     \-- ...
    def get_connectivity(self, interactions):
        _connectivity = dict()
        for key in interactions.keys():
            if interactions[key]['N'] == 0:
                self.logger.warn(f'No {key} information found in GSD file.')
                _connectivity[key] = None
            else:
                _connectivity[key] = list(
                    map(tuple, interactions[key]['group'].tolist())
                )

        _connectivity['particles_group'] = self.get_particles_group(_connectivity)

        return _connectivity

    def get_system_info(self, frame_idx=None, frame=None):
        self._system_info = {'system': dict(), 'outputs': dict()}
        _path = f'{frame_idx}'
        _interaction_types = [
            'bonds',
            'angles',
            'dihedrals',
            'impropers',
            'constraints',
            'pairs',
        ]

        def get_particle_parameters(path: str = None, frame=None):
            n_particles = self._data_parser.get(f'{path}.N', frame=frame)
            if n_particles is None:
                return dict()
            else:
                return self._data_parser.get(path, frame=frame).__dict__

        self._particle_data_dict = get_particle_parameters(
            path=f'{frame_idx}.particles', frame=frame
        )

        if self._particle_data_dict is None:
            self.logger.warning(
                f'No number of particles available in frame {frame_idx}. Other'
                ' particle attributes will not be stored for frame {frame_idx}.'
            )

        def get_value(value, path=None, frame=None):
            if value is None:
                return value
            try:
                value = self._data_parser.get(
                    f'{path}.{value}' if path else value, frame=frame
                ).__dict__
                if value is None:
                    self.logger.warning(
                        f'No attributes found for key {value}. {value.upper()} attributes will '
                        'not be stored.'
                    )
                    return None
                else:
                    return value
            except KeyError:
                return value
            # TODO: handle data chunks that are not standard objects
            except AttributeError:
                pass

        info_keys = {
            'step': ['system', 'outputs'],
            'n_atoms': 'system',
            'positions': 'system',
            'labels': 'system',
            'mass': 'system',
            'velocities': 'system',
            'charge': 'system',
            'lattice_vectors': 'system',
            'periodic_boundary_conditions': 'system',
            'dimensionality': 'system',
            'forces': 'outputs',
            'label': 'outputs',
        }

        # Get quantities from particles chunk of GSD file
        for key, gsd_key in self._nomad_to_particles_group_map.items():
            section = info_keys[key]
            if isinstance(section, list):
                for sec in section:
                    self._system_info[sec][key] = (
                        self._particle_data_dict[gsd_key]
                        if gsd_key is not None
                        else None
                    )
            else:
                self._system_info[section][key] = (
                    self._particle_data_dict[gsd_key] if gsd_key is not None else None
                )
        # Get step and box attributes from configurations chunk of GSD file:
        for key, gsd_key in self._nomad_to_box_group_map.items():
            section = info_keys[key]
            _values_dict = get_value('configuration', path=_path, frame=frame)
            if isinstance(section, list):
                for sec in section:
                    self._system_info[sec][key] = (
                        _values_dict[gsd_key] if gsd_key is not None else None
                    )
            else:
                self._system_info[section][key] = (
                    _values_dict[gsd_key] if gsd_key is not None else None
                )

        # Extract interacton infromation from frame,
        # build connectivity structure following Nomad-H5MD schema.
        _interaction_dicts = dict()
        for interaction in _interaction_types:
            _interaction_dicts[interaction] = get_value(
                interaction, path=_path, frame=frame
            )
        self._connectivity = self.get_connectivity(_interaction_dicts)
        # print(self._connectivity)

        return self._system_info

    def parse_system(self, simulation, frame_idx=None, frame=None):
        system_info = self._system_info.get(
            'system'
        )  # ? Will I really need this in the end?
        atoms_dict = self._system_info.get('system')
        _path = f'{frame_idx}'
        if not system_info:
            self.logger.error('No particle information found in GSD file.')
            return

        self._system_time_map = {}  # ? Is this required?
        topology = None  # ? Is this required?

        # ! This was added according to existing nomad approach,
        # ! does not support (semi) grand canonical hoomd-Blue output.
        if self._first_frame is True:
            self.logger.warn(
                'Only the topology of the first frame will be stored, '
                'grand canonical simulations are currently not supported.'
            )
            atoms_dict['is_representative'] = True
            self._first_frame = False
        else:
            atoms_dict['is_representative'] = False

        atom_labels = atoms_dict.get('labels')
        if atom_labels is not None:
            try:
                symbols2numbers(atom_labels)
                atoms_dict['labels'] = atom_labels
            except KeyError:  # TODO this check should be moved to the system normalizer in the new schema
                atoms_dict['labels'] = ['X'] * len(atom_labels)

        bond_dict = self._data_parser.get(f'{frame_idx}.bonds', frame=frame).__dict__
        atoms_dict['bond_list'] = bond_dict['group']

        # ! Natively, no time step stored in GSD file. Copy frame index instead,
        # ! alert user to missing information.
        time = atoms_dict.pop('step')
        time_unit = time.units if hasattr(time, 'units') else None
        atoms_dict['time_step'] = time.magnitude if time_unit is not None else time
        if time_unit is None:
            self.logger.warn(
                'No magnitude and unit information provided for the '
                'simulation time step'
            )

        # path_topology = 'connectivity.particles_group'
        # topology = self._data_parser.get(path_topology)

        # # REMAP some of the data for the schema
        # atoms_dict['branch_label'] = (
        #     'Total System'  # ? Do we or should we have a default name for the entire system?
        # )
        # atoms_dict['time_step'] = atoms_dict.pop(
        #     'time'
        # ).magnitude  # TODO change in system_info
        # atomic_cell_keys = [
        #     'n_atoms',
        #     'lattice_vectors',
        #     'periodic_boundary_conditions',
        #     'positions',
        #     'velocities',
        #     'labels',
        # ]
        # atoms_dict['atomic_cell'] = {}
        # for key in atomic_cell_keys:
        #     atoms_dict['atomic_cell'][key] = atoms_dict.pop(key)

        # self.parse_trajectory_step(atoms_dict, simulation)

        # if i_step == 0 and topology:  # TODO extend to time-dependent topologies
        #     self.parse_system_hierarchy(
        #         simulation.model_system[-1], topology, path_topology
        #     )

    def write_to_archive(self) -> None:
        #######################################################################
        # Access simulation file(s).
        #######################################################################
        self._maindir = os.path.dirname(
            self.mainfile
        )  # ? GSD output single file or more?
        self._gsd_files = [
            _file for _file in os.listdir(self._maindir) if _file.endswith('.gsd')
        ]
        self._basename = os.path.basename(self.mainfile).rsplit('.', 1)[0]
        self._data_parser.mainfile = self.mainfile
        if self._data_parser.filegsd is None:
            self.logger.warning('GSD file missing in GSD Parser.')
            return

        #######################################################################
        # Start populating NOMAD-simulations schema
        #######################################################################
        simulation = Simulation()
        self._program_dict = self.get_program_info()
        simulation.program = Program(
            name=self._program_dict['name'],
            version=self._program_dict['version'],
        )
        for frame_idx, frame in enumerate(self._data_parser.filegsd):
            self.get_system_info(frame_idx=frame_idx, frame=frame)
            self.parse_system(simulation, frame_idx=frame_idx, frame=frame)
            # print(self._system_info)
            # TODO: parse systems_info dict to nomad

            # ? Forces etc. are user-defined and read via get_logged_info?
            # TODO: Extract observables from logged data, parse to ModelOutput
            # info_keys = {
            #     'forces': 'calculation',
            # }
