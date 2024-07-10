import sys

def process_scaffold(fastafile, scaffoldname, output):
    found_scaffold = False
    with open(fastafile, 'r') as infile:
        with open(output, 'w') as outfile:
            current_sequence = []
            for line in infile:
                if line.startswith('>'):
                    # If we find a header line
                    if scaffoldname in line:
                        found_scaffold = True
                    else:
                        found_scaffold = False
                        continue
                if found_scaffold:
                    outfile.write(line)

# Main function to handle command-line arguments
def main():
    if len(sys.argv) != 4:
        print("Usage: python3 ScaffoldExtractor.py fastafile scaffoldname output")
        return

    fastafile = sys.argv[1]
    scaffoldname = sys.argv[2]
    output = sys.argv[3]
    
    process_scaffold(fastafile, scaffoldname, output)

if __name__ == "__main__":
    main()

