from pathlib import Path

# files = input("Enter file name(s), separated by commas: ")
# input_file_paths = files.split(",")
# print(f"You put in {len(input_file_paths)} files.")
count = 1
# What accession would you like to align? 
accessions = ["PC185"]

# The next line sets the genome you want to align that accession to 
input_file_paths = ["L58_LIBRARIES.txt"]
for mappingTo in accessions:    
    # Define the input file path
    for file in input_file_paths:
        input_file_path = "lib/"+file
        file_name = file[:-14]
        to_clean = "rawCounts/"+file_name+"_tpm_counts.txt"
        outputDir = "out/" + mappingTo + "/"
        Path(outputDir).mkdir(parents=True, exist_ok=True)
        output = outputDir + file_name+"_BestAlignment_tpm_counts" + mappingTo + ".txt"

        # Initialize a list to store the gene library names that do not meet the criteria
        gene_library_names = []
        names_to_change = {}

        # Open the input file for reading
        with open(input_file_path, "r") as input_file:
            # Skip the header line
            next(input_file)

            # Iterate through the dataset and filter rows
            for line in input_file:
                columns = line.strip().split('\t')
            
                groupName = columns[7]  # Column 8 (groupName)
                sampleName = columns[4]  # Column 5 (sampleName)

                # Check if the groupName is "filtered reads" or sampleName contains "A03"
                if groupName == "filtered reads" or mappingTo not in sampleName:
                    gene_library_names.append(columns[0])  # Append the gene library name (Column 1)
                else:
                    names_to_change[columns[0]] = columns[4]

        # Open the input file for reading
        with open(to_clean, "r") as input_file:
            # Read the header line and split it into a list of column headers
            headers = input_file.readline().strip().split('\t')  # Assuming tab-separated columns

            # Determine the indices of columns to keep
            indices_to_keep = [i for i, header in enumerate(headers) if header not in gene_library_names]
            indices_to_change = [i for i, header in enumerate(headers) if header in names_to_change]

            # Create a new list of headers with only the columns to keep
            new_headers = [headers[i] for i in indices_to_keep]
            for i, header in enumerate(new_headers):
                if header in names_to_change:
                    new_headers[i] = names_to_change[header]

            # Create a new file for writing the modified dataset
            with open(output, "w") as output_file:
                # Write the modified header line to the output file
                output_file.write('\t'.join(new_headers) + '\n')

                # Process and write the remaining lines
                for line in input_file:
                    values = line.strip().split('\t')
                    new_values = [values[i] for i in indices_to_keep]
                    output_file.write('\t'.join(new_values) + '\n')
        print(f"File '{output}' is done, {len(input_file_paths)-count} files left!")
        count+=1
    count = 1

