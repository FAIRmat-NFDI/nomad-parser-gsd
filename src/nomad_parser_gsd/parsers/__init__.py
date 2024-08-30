from nomad.config.models.plugins import ParserEntryPoint
from pydantic import Field


class GSDParserEntryPoint(ParserEntryPoint):
    parameter: int = Field(0, description='Custom configuration parameter')

    def load(self):
        from nomad_parser_gsd.parsers.parser import GSDParser

        return GSDParser(**self.dict())


parser_entry_point = GSDParserEntryPoint(
    name='GSDParser',
    description='Parser for the GSD file format, the native file format for HOOMD-blue. (https://gsd.readthedocs.io/en/v3.3.1/, https://glotzerlab.engin.umich.edu/hoomd-blue/)',
    mainfile_name_re='.*\.gsd.*',
    mainfile_contents_re=r'GSD version [\d\.]*',
)
