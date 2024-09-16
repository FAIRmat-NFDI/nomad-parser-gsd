import os
from typing import TYPE_CHECKING, Iterable, Union

import numpy as np

# from nomad.parsing.parser import MatchingParser
import structlog
from atomisticparsers.utils import MDParser
from nomad.config import config
from nomad.datamodel.metainfo.workflow import Workflow
from nomad.parsing.file_parser import FileParser
from nomad_simulations.schema_packages.general import (
    Program,
    Simulation,
)
from nomad.units import ureg

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
            except Exception:
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
        self._time_unit = ureg.picosecond

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

        return _program_info

    # Additional data are stored in the log dictionary as numpy arrays:
    """
    Logged data encompasses values computed at simulation time that are too expensive 
    or cumbersome to re-compute in post processing. This specification does not define 
    specific chunk names or define logged data. Users may select any valid name for 
    logged data chunks as appropriate for their workflow.
    """

    def get_logged_info(self):
        try:
            return gsdhoomd.read_log(name=self.mainfile, scalar_only=False)
        except FileNotFoundError:
            self.logger.warning(
                'No additional logged data found, no user-defined data will be stored.'
            )
            return dict()

    # def get_configuration_info(self, path=None, frame=None):
    #     if not frame:  # TODO: come up with useful sanity check
    #         return dict()
    #     else:
    #         return self._data_parser.get(path, frame=frame).__dict__

    def get_particle_parameters(self, path: str = None, frame=None):
        n_particles = self._data_parser.get(f'{path}.N', frame=frame)
        if n_particles is None:
            return dict()
        else:
            return self._data_parser.get(path, frame=frame).__dict__

    def get_system_info(self):
        self._system_info = {'system': dict(), 'calculation': dict()}
        self._n_particles_all = np.array([], dtype=int)
        self._particles_positions_all = list()
        self._particles_types_all = list()
        self._velocities_all = list()
        self._simulation_steps_all = np.array([], dtype=int)
        self._dimensions_all = np.array([], dtype=int)
        self.box_all = list()  # np.array((self.n_frames, 6), dtype=float)
        n_frames = self._n_frames

        # def get_value(value, steps, path=None):
        #     if value is None:
        #         return value
        #     try:
        #         value = self._data_parser.get(f'{path}.value' if path else 'value')
        #         if value is None:
        #             self.logger.warning(
        #                 'Missing values in particle attributes.'
        #                 ' These attributes will not be stored.'
        #             )
        #             return None
        #         else:
        #             return value
        #     except KeyError:
        #         return [
        #             value
        #         ] * n_frames  # TODO: what default makes sense if key not populated?

        for frame_idx, frame in enumerate(self._data_parser.filegsd):
            _particle_data_dict = self.get_particle_parameters(
                path=f'{frame_idx}.particles', frame=frame
            )
            if _particle_data_dict is None:
                self.logger.warning(
                    f'No number of particles available in frame {frame_idx}. Other'
                    ' particle attributes will not be stored for frame {frame_idx}.'
                )
            self._n_particles_all = np.append(
                self._n_particles_all, _particle_data_dict['N']
            )
            self._particles_positions_all.append(_particle_data_dict['position'])
            self._particles_types_all.append(_particle_data_dict['types'])
            self._velocities_all.append(_particle_data_dict['velocity'])

            _configuration_data_dict = self._data_parser.get(
                f'{frame_idx}.configuration', frame=frame
            ).__dict__

            self._simulation_steps_all = np.append(
                self._simulation_steps_all, _configuration_data_dict['step']
            )
            _bonds_dict = self._data_parser.get(
                f'{frame_idx}.bonds', frame=frame
            ).__dict__
            _angles_dict = self._data_parser.get(
                f'{frame_idx}.angles', frame=frame
            ).__dict__
            _dihedrals_dict = self._data_parser.get(
                f'{frame_idx}.dihedrals', frame=frame
            ).__dict__
            _impropers_dict = self._data_parser.get(
                f'{frame_idx}.impropers', frame=frame
            ).__dict__
            _special_pairs_dict = self._data_parser.get(
                f'{frame_idx}.pairs', frame=frame
            ).__dict__

        self._system_info['system']['positions'] = np.array(
            self._particles_positions_all
        )
        self._system_info['system']['n_atoms'] = self._n_particles_all
        self._system_info['system']['labels'] = np.array(self._particles_types_all)
        self._system_info['system']['velocities'] = np.array(self._velocities_all)

        self._system_info['calculation']['steps'] = self._simulation_steps_all

        return self._system_info

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

        self.system_info = self.get_system_info()
        print(self.system_info)
        # TODO: parse systems_info dict to nomad

        # ? Forces etc. are user-defined and read via get_logged_info?
        # system_keys = {
        #     'forces': 'calculation',
        # }
