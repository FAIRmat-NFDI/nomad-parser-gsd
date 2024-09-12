from typing import TYPE_CHECKING, Union, Iterable

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )
    from structlog.stdlib import (
        BoundLogger,
    )

from atomisticparsers.utils import MOL, MDParser
from nomad.config import config
import os
from nomad.datamodel.metainfo.workflow import Workflow
from nomad.parsing.file_parser import FileParser
import structlog
# from nomad.parsing.parser import MatchingParser

logging = structlog.get_logger()
try:
    import gsd.hoomd as gsdhoomd
    # import gsd.pygsd as gsdpy
except ImportError:
    logging.warn('Required module gsd.hoomd not found.')
    gsdhoomd = False
    # gsdpy = False

configuration = config.get_plugin_entry_point(
    'nomad_parser_gsd.parsers:parser_entry_point'
)

MOL = 6.022140857e23


class GSDFileParser(FileParser):  # or MatchingParser?
    def __init__(self):  # , *args, **kwargs
        super().__init__(None)  # , *args, **kwargs

    #         self._nomad_to_hoomdblue_map = {}
    #         self._nomad_to_hoomdblue_map['system'] = {}
    #         self._nomad_to_hoomdblue_map['system']['atoms'] = {
    #             # 'lattice_vectors': 'configuration.box',
    #             'positions': 'particles.position',
    #             'x_hoomdblue_orientation': 'particles.orientation',
    #             # 'x_hoomdblue_typeid': 'particles.typeid',
    #             'velocities': 'particles.velocity',
    #             'x_hoomdblue_angmom': 'particles.angmom',
    #             'x_hoomdblue_image': 'particles.image',
    #         }

    #         self._nomad_to_hoomdblue_map['calculation'] = {
    #             'step': 'configuration.step',
    #         }

    #         self._nomad_to_hoomdblue_map['method'] = {}
    #         self._nomad_to_hoomdblue_map['method']['atom_parameters'] = {atomistic-parsers
    #             # 'x_hoomdblue_types': 'particles.types',
    #             'mass': 'particles.mass',
    #             'charge': 'particles.charge',
    #             'x_hoomdblue_diameter': 'particles.diameter',
    #             'x_hoomdblue_body': 'particles.body',
    #             'x_hoomdblue_moment_inertia': 'particles.moment_inertia',
    #             'x_hoomdblue_type_shapes': 'particles.type_shapes',
    #         }
    #         self._hoomdblue_interaction_keys = [
    #             'bonds',
    #             'angles',
    #             'dihedrals',
    #             'impropers',
    #             'constraints',
    #             'pairs',
    #         ]

    #     def get_attribute(self, source, path: str = None, default=None):
    #         """
    #         Extracts attribute from object based on path, and returns default if not defined.
    #         """
    #         if path:
    #             section_segments = path.split('.')
    #             for section in section_segments:
    #                 try:
    #                     value = getattr(source, section)
    #                     source = value[-1] if isinstance(value, list) else value
    #                 except Exception:
    #                     return
    #             source = source if source is not None else default
    #             return source

    #     def apply_unit(self, quantity, unit: str, unit_factor: float):
    #         if quantity is None:
    #             return
    #         if unit:
    #             unit_val = ureg(unit)
    #             unit_val *= unit_factor
    #             quantity *= unit_val

    #         return quantity

    def filegsd(self):
        print(self._file_handler)
        if self._file_handler is None:
            print(self.mainfile)
            try:
                self._file_handler = gsdhoomd.open(name=self.mainfile, mode='r')
            except Exception:
                self.logger.error('Error reading gsd file.')
        return self._file_handler

    def parse(self, path: str = None, **kwargs):
        source = kwargs.get('source', self.filegsd)
        isattr = kwargs.get('isattr', False)
        print(source, isattr)
        # value = None
        # if isattr:
        #     attr_path, attribute = path.rsplit('.', 1)
        #     value = self.get_attribute(source, attribute, path=attr_path)
        # # else:
        # #     value = self.get_value(source, path)
        # self._results[path] = value


class GSDParser(MDParser):
    def __init__(self, *args, **kwargs):
        super().__init__()
        self._data_parser = GSDFileParser()
        self._maindir = None
        self._gsd_files = None
        self._basename = None

    def test(self):
        for frame in enumerate(self._data_parser.filegsd):
            print(frame)

    def write_to_archive(self) -> None:
        self._maindir = os.path.dirname(
            self.mainfile
        )  # GSD output single file or more?
        self._gsd_files = [
            _file for _file in os.listdir(self._maindir) if _file.endswith('.gsd')
        ]
        self._basename = os.path.basename(self.mainfile).rsplit('.', 1)[0]
        self._data_parser.mainfile = self.mainfile

        print('_maindir:', self._maindir)
        print('_gsd_files:', self._gsd_files)
        print('_basename:', self._basename)
