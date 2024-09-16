from typing import TYPE_CHECKING, Union, Iterable

import numpy as np

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
from nomad.units import ureg
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

    # Extract data from gsd file, store with keyword (to_nomad mapping here or downstairs?)
    # Keep sections separate -> easier to handle/debug

    def filegsd(self):
        if self._file_handler is None:
            try:
                self._file_handler = gsdhoomd.open(name=self.mainfile, mode='r')
            except Exception:
                self.logger.error('Error reading gsd file.')
        return self._file_handler

    def get_attribute(self, source, path: str = None, default=None):
        """
        Extracts attribute from object based on path, and returns default if not defined.
        """
        if path:
            section_segments = path.split('.')
            print(section_segments)
            for section in section_segments:
                print(section)
                try:
                    value = getattr(source, section)
                    print(value)
                    source = value[-1] if isinstance(value, list) else value
                except Exception as e:
                    print(e)
                    return
            source = source if source is not None else default
            print(source)
            return source

    def parse(self, path: str = None, **kwargs):
        pass
        # print('kwargs', kwargs)
        # source = kwargs.get('source', self.filegsd())
        # isattr = kwargs.get('isattr', False)
        # print(source, isattr)
        # value = None
        # if isattr:
        #     attr_path, attribute = path.rsplit('.', 1)
        #     print(attr_path, attribute)
        #     value = self.get_attribute(source, attribute, path=attr_path)
        #     print(value)
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
        self._n_frames = None
        self._n_atoms = None
        self._atom_parameters = None
        self._time_unit = ureg.picosecond
        self._frame_configuration_step = 'configuration.step'
        self._frame_configuration_dimensions = 'configuration.dimensions'
        self._frame_configuration_box = 'configuration.box'
        self._frame_particles_position = 'particles.position'

    # Populating nomad schema and writing to archive here
    def get_system_info(self, frame) -> dict:
        self._system_info = {'system': {}, 'calculation': {}}

    def write_to_archive(self) -> None:
        self._maindir = os.path.dirname(
            self.mainfile
        )  # GSD output single file or more?
        self._gsd_files = [
            _file for _file in os.listdir(self._maindir) if _file.endswith('.gsd')
        ]
        self._basename = os.path.basename(self.mainfile).rsplit('.', 1)[0]
        self._data_parser.mainfile = self.mainfile
        if self._data_parser.filegsd is None:
            self.logger.warning('GSD file missing in GSD Parser.')
            return

        for frame_idx, frame in enumerate(self._data_parser.filegsd()):
            print(len(self._data_parser.filegsd()))
            print(frame_idx, frame)
            configuration = self._data_parser.get_attribute(frame, 'configuration')
            print('configuration', configuration)
            # positions = self._data_parser.get(frame, 'particles.positions', None)
            # print('positions', positions)
