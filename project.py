# This file reads in a CSV and parses it to get a list of all RSIDs represnted as edges,
# convert them to genes with the MyVariant API, and then re-write the CSV with the genes.

import csv
import myvariant

# Create a MyVariant API object
mv = myvariant.MyVariantInfo()

# Make a (global) result map between RSID and gene(s)
resultMap = {}

inputPath = input("Enter the name of the input graph file: ")
outputPath = input("Enter the desired name of the output graph file: ")

# Read in the CSV
with open(inputPath, 'r+') as f:
    csvReader = csv.DictReader(f)

    # Get all RSIDs
    RSIDs = set()
    for line in csvReader:
        for rsID in line['Label'].split(","):
            RSIDs.add(rsID)
    
    # Convert RSIDs to genes
    results = mv.querymany(list(RSIDs), scopes='dbsnp.rsid')

    for result in results:
        # Check if there are any errors
        if not isinstance(result, dict):
            print("Result is not a dict: ")
            print(result)
            continue
        if 'error' in result:
            print("Error:")
            print(result['error'])
            continue
        if 'notfound' in result:
            print("Not found:")
            print(result['notfound'])
            continue
        if 'dbsnp' not in result:
            print("No dbsnp:")
            print(result)
            continue
        if 'query' not in result:
            print("No query:")
            print(result)
            continue
        if 'gene' not in result['dbsnp']:
            print("Intergenic SNP", result['query'], "found!")
            # Means that it is an intergenic SNP, so will have multiple genes (possibly)
            # The parenthesis show that it is an intergenic SNP
            resultMap[result['query']] = ",".join(result['snpeff']['ann']['genename'].split("-"))
        else:
            # Check if the result['dbsnp']['gene'] is a list or a dict
            if isinstance(result['dbsnp']['gene'], list):
                for gene in result['dbsnp']['gene']:
                    if result["query"] in resultMap:
                        resultMap[result['query']][gene['symbol']] = gene['geneid']
                    else:
                        resultMap[result['query']] = {gene['symbol'] : gene['geneid']}
            else:
                gene = result['dbsnp']['gene']
                if result["query"] in resultMap:
                    resultMap[result['query']][gene['symbol']] = gene['geneid']
                else:
                    resultMap[result['query']] = {gene['symbol'] : gene['geneid']}

    for rsID, result in resultMap.items():
        print(rsID, ":", result)

    # Re-write the CSV with the genes
    f.seek(0)
    csvReader = csv.DictReader(f)

    # Create a new CSV file
    with open(outputPath, 'w+') as newCSV:
        fieldnames = ['Source', 'Target', 'Type', 'Id', 'Label', 'Weight']
        csvWriter = csv.DictWriter(newCSV, fieldnames=fieldnames)
        csvWriter.writeheader()

        # Write the new CSV
        for line in csvReader:
            # Collect all genes for this edge
            genes = []
            seenGenes = set()
            for rsID in line['Label'].split(","):
                # Check if it is a dict (meaning on gene) or a list (meaning intergenic)
                if isinstance(resultMap[rsID], dict):
                    for gene, geneID in resultMap[rsID].items():
                        if gene not in seenGenes:
                            seenGenes.add(gene)
                            genes.append(gene)
                elif resultMap[rsID] not in seenGenes:
                    seenGenes.add(resultMap[rsID])
                    genes.append("(" + resultMap[rsID] + ")")
            # Now write the new line
            csvWriter.writerow({
                'Source': line['Source'],
                'Target': line['Target'],
                'Type': line['Type'],
                'Id': line['Id'],
                'Label': ",".join(genes),
                'Weight': len(genes)
            })
        # Close the new CSV
        newCSV.close()
    # Close the old CSV
    f.close()

print("Done!")
