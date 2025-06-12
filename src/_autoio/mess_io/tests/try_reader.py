import sys
import mess_io.reader.inp as inp
import ioformat.pathtools as parser
from ioformat import remove_comment_lines
import mess_io.reader._mess_keys as KEYS
import autoparse.pattern as app
import autoparse.find as apf
from autoparse import cast as ap_cast

#pattern = app.capturing(app.NUMBER) + app.one_of_these([app.SPACES, app.LINE_END])
#print(pattern)
#test_str = 'C 1 2 3 4'
#captures = apf.all_captures(pattern, test_str)
#print(captures)
#print(list(ap_cast(captures)))

mess_filename = 'data/ME2s_a.inp'

JOB_PATH = sys.argv[1]

input_str = parser.read_file(JOB_PATH, mess_filename)
input_lines = input_str.split('\n')

model_lines, matches, global_lines = inp.get_sections(input_lines, KEYS.OVERALL)
model_sections, matches, _ = inp.get_sections(model_lines[0], KEYS.MODEL)

# Parse through each section
for idx, model_section in enumerate(model_sections):
    print('Model section match:\n', matches[idx])
    #print('section:\n', model_section)
    inp.parse_model_section(model_section, matches[idx])


