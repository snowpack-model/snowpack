###############################################################################
#   Copyright 2022 WSL Institute for Snow and Avalanche Research  SLF-DAVOS   #
###############################################################################
# This file is part of INIshell.
# INIshell is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# INIshell is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with INIshell.  If not, see <http://www.gnu.org/licenses/>.


"""
Validation of XML files based on a given XSD schema

Usage:
Pass the path to the XSD followed by an arbitrary amount of paths to XML files

Example: 
$ python3 validate.py schema.xsd file.xml file2.xml
"""

import sys
from typing import List
from xsd_validator import XsdValidator, XsdValidationError


class Validator:

    def __init__(self, xsd_file: str):
        self.xsd_validator = XsdValidator(xsd_file)

    def validateXML(self, xml_file: str) -> bool:
        try:
            print("Validating {0}".format(xml_file))
            self.xsd_validator.assert_valid(xml_file)
            print("Validation successful")
            return True

        except XsdValidationError as err:
            print("Validation failed: {0}".format(err))
            return False

        except OSError as err:
            print(err)
            return False


def main(xsd_file: str, xml_files: List[str]) -> bool:
    validator = Validator(xsd_file=xsd_file)
    errors = 0
    for xml_file in xml_files:
        errors += 0 if validator.validateXML(xml_file=xml_file) else 1
    print("Validated {0} files and found {1} error(s)".format(
        len(xml_files), errors))
    exit(0 if errors == 0 else 1)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Please provide a path to an XSD file followed by an arbitrary amount of paths to XML files")
        print("Example: $ python3 validate.py schema.xsd file.xml file2.xml")
        exit(1)
    main(sys.argv[1], sys.argv[2:])
