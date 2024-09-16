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


class GSDFileParser(FileParser):  # or MatchingParser?
    def __init__(self):  # , *args, **kwargs
        super().__init__(None)  # , *args, **kwargs

    # Extract data from gsd file, store with keyword (to_nomad mapping here or downstairs?)
    # Keep sections separate -> easier to handle/debug

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

    # def get_attribute(self, source, attribute: str = None, default=None):
    #     """
    #     Extracts attribute from object based on path, and returns default if not defined.
    #     """
    #     if attribute:
    #         section_segments = attribute.split('.')
    #         for section in section_segments:
    #             try:
    #                 value = getattr(source, section)
    #                 source = value[-1] if isinstance(value, list) else value
    #             except Exception:
    #                 return
    #         source = source if source is not None else default
    #         return source

    def get_value(self, group, path: str, default=None):
        """
        Extracts group or dataset from group object based on path, and returns default if not defined.
        """
        section_segments = path.split('.')
        for section in section_segments:
            try:
                value = getattr(group, section)  # group.get(section)
                # unit = self.get_attribute(group, 'unit', path=section)
                # unit_factor = self.get_attribute(
                #     group, 'unit_factor', path=section, default=1.0
                # )
                group = value
            except Exception:
                return

        # if value is None:
        #     value = default
        # elif isinstance(value, h5py.Dataset):
        #     value = value[()]
        #     value = self.apply_unit(value, unit, unit_factor)
        # value = self.decode_bytes(value)

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
        self._n_frames = None
        self._n_atoms = None
        self._atom_parameters = None
        self._time_unit = ureg.picosecond
        self._frame_particles_position = 'frame.particles.position'
        # self._particles_position_all = np.array([])

    # Populating nomad schema and writing to archive here

    def write_to_archive(self) -> None:
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

        for frame_idx, frame in enumerate(self._data_parser.filegsd):
            positions = self._data_parser.get(
                f'{frame_idx}.particles.position', frame=frame
            )
            print(positions)
