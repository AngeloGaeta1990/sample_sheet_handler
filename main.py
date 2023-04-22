import sys
import argh
from sample_sheet_handler import SampleSheet
from sample_sheet_handler import reverse_complement_index2
from sample_sheet_handler import convert_sample_sheet_to_v2




parser = argh.ArghParser()
parser.add_commands([reverse_complement_index2, convert_sample_sheet_to_v2])

if __name__ == "__main__":
    parser.dispatch()