import argparse
parser = argparse.ArgumentParser()
parser.add_argument("path", help="The path to the input CSV file")
parser.add_argument("-d","--debug", help="Runs the code in debug mode", action="store_true")
parser.add_argument("-g","--graphs", help="Outputs the graph", action="store_true")
args = parser.parse_args()

print(args.path)

if args.debug:
	print("debug mode")
else:
	print("normal mode")
	
	
if args.graphs:
	print("----PICTURE----")