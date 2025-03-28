
import sys
import re
import random
import matplotlib.pyplot as plt

## USAGE ##
# python3 multiple_alignments.py chr[n]

# Parse GTF file to get coding regions
def parse_gtf(gtf_ifn):
    coding_regions = []
    with open(gtf_ifn, "r") as gtf:

        last_field = ''
        strand = ''
        last_cds_search = False

        for line in gtf:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if fields[2] == 'transcript':
                if last_cds_search:
                    start, end, strand, gene_name = temp
                    rc_start = end
                    rc_end = start
                    coding_regions.append((rc_start, rc_end, strand, gene_name))
                    last_cds_search = False

                strand = fields[6]
            if strand == '+':
                if fields[2] == "CDS" and (last_field == 'transcript' or last_field == 'exon'):
                    start, end = int(fields[3]), int(fields[4])
                    items = fields[8].split(";")
                    gene_name = items[2].strip()
                    coding_regions.append((start, end, strand, gene_name))
            elif strand == '-':
                if fields[2] == "CDS":
                    last_cds_search = True
                    start, end = int(fields[3]), int(fields[4])
                    items = fields[8].split(";")
                    gene_name = items[2].strip()
                    temp = (start, end, strand, gene_name)

            last_field = fields[2]
    return coding_regions

# Expand coding regions to left and right of translation initiation site
def get_expanded_regions(coding_regions):

    expanded_coding_regions = []
    flank = 500
    for (start, end, strand, gene_name) in coding_regions:
        if strand == '+':
            expanded_start = start - flank
            expanded_end = start + flank
            expanded_coding_regions.append([expanded_start, expanded_end, start, end, strand, gene_name])
        elif strand == '-':
            expanded_start = start + flank
            expanded_end = start - flank
            expanded_coding_regions.append([expanded_start, expanded_end, start, end, strand, gene_name])

    return expanded_coding_regions

def get_neg_limit(chrom):
    if chrom == 'chr1':
        n_genes = 3000
    elif chrom == 'chr2':
        n_genes = 2500
    elif chrom == 'chr3':
        n_genes = 1900
    elif chrom == 'chr4':
        n_genes = 1600
    elif chrom == 'chr5':
        n_genes = 1700
    elif chrom == 'chr6':
        n_genes = 1900
    elif chrom == 'chr7':
        n_genes = 1800
    elif chrom == 'chr8':
        n_genes = 1400
    elif chrom == 'chr9':
        n_genes = 1400
    elif chrom == 'chr10':
        n_genes = 1400
    elif chrom == 'chr11':
        n_genes = 2000
    elif chrom == 'chr12':
        n_genes = 1600
    elif chrom == 'chr13':
        n_genes = 800
    elif chrom == 'chr14':
        n_genes = 1200
    elif chrom == 'chr15':
        n_genes = 1200
    elif chrom == 'chr16':
        n_genes = 1300
    elif chrom == 'chr17':
        n_genes = 1600
    elif chrom == 'chr18':
        n_genes = 600
    elif chrom == 'chr19':
        n_genes = 1700
    elif chrom == 'chr20':
        n_genes = 900
    elif chrom == 'chr21':
        n_genes = 400
    elif chrom == 'chr22':
        n_genes = 800 
    return n_genes * 4

def get_pos_sites(expanded_coding_regions, start_pos, end_pos):
    in_bounds = False
    cds_idx = -1
    
    for (expanded_start, expanded_end, true_start, true_end, ref_strand, gene_name) in expanded_coding_regions:
        if ref_strand == '+':
            if ((start_pos >= expanded_start and start_pos < expanded_end) 
                or (start_pos < expanded_start and end_pos > expanded_start) 
                or (start_pos < expanded_end and end_pos > expanded_end)):
                in_bounds = True
                if true_start >= start_pos and true_start < end_pos:
                    cds_idx = true_start - start_pos - 1
                break
        elif ref_strand == '-':
            if ((start_pos >= expanded_end and start_pos < expanded_start) 
                or (start_pos < expanded_end and  end_pos > expanded_end) 
                or (start_pos < expanded_start and end_pos > expanded_start)):
                in_bounds = True
                if true_start >= end_pos and true_start < start_pos:
                    cds_idx = start_pos - true_start - 1
                break
    
    return in_bounds, cds_idx, gene_name

# Create dictionary to store gene sequences
def create_dicts(maf_ifn, expanded_coding_regions, chrom):
    pos_dict = {}
    neg_dict = {}
    pos_idx = 0
    neg_idx = 0
    neg_count = 0
    neg_gene_name = ''

    set_limit = get_neg_limit(chrom)
    
    with open(maf_ifn, "r") as maf:
        for line in maf:

            if line.startswith("a"):  # Alignment block header
                subseqs = dict()
                in_bounds = True
                skip_neg = False

            elif line.startswith("s"):
                if skip_neg: continue

                fields = re.split(r"\s+", line)
                chrom = fields[1]
                strand = fields[4]
                seq = fields[6]

                # Find start and end positions
                if "hg38" in line:                  # adjust for human genome identifier

                    start = int(fields[2])
                    size = int(fields[3])
                    src_size = int(fields[5])

                    start_pos = start if strand == "+" else start + size
                    end_pos = start + size if strand == "+" else start

                    in_bounds, cds_idx, gene_name = get_pos_sites(expanded_coding_regions, start_pos, end_pos)
                   
                    # Stop collecting negative control examples at limit
                    if not in_bounds and neg_count > set_limit and random.randint(1, 100) > 1:
                        skip_neg = True

                    # Update negative control idx and gene name
                    if cds_idx != -1:     # If in CDS segment, skip to next negative control gene name
                        neg_count += 1
                        while (neg_idx%10) != 0:
                            neg_idx+=1
                    elif in_bounds and neg_count > set_limit and (neg_idx%10) == 0:
                        neg_count += 1
                    else:
                        if neg_count <= set_limit and (neg_idx%10) == 0:
                            neg_count += 1

                    neg_gene_name = f"NEG_CTRL_{neg_count}"
                        

                subseqs[chrom] = [strand, seq]

            elif line.strip() == "":  # End of block
                if cds_idx != -1:   # in cds segment
                    pos_dict[pos_idx] = [start_pos, end_pos, size, gene_name, cds_idx, subseqs]
                    pos_idx += 1
                elif in_bounds:     # in bounds but not in cds segment
                    pos_dict[pos_idx] = [start_pos, end_pos, size, gene_name, cds_idx, subseqs]
                    pos_idx += 1
                    if neg_count > set_limit:
                        neg_dict[neg_idx] = [start_pos, end_pos, size, neg_gene_name, cds_idx, subseqs]
                        neg_idx += 1
                else:               # not in bounds, only negative
                    if skip_neg: continue
                    neg_dict[neg_idx] = [start_pos, end_pos, size, neg_gene_name, cds_idx, subseqs]
                    neg_idx += 1        

    print("create_dicts done")
    return pos_dict, neg_dict

def process_pos_dict(pos_dict, expanded_coding_regions, ofn):
    
    with open(ofn, "w") as out:
        
        gene_dict = dict()
        # {gene_a: [leftmost, rightmost, cds_block_start, cds_idx, {human: actg, chimp: actt, ...}] }

        curr_gene = ''
        for idx in pos_dict:
            start_pos, end_pos, size, gene_name, cds_idx, subseqs = pos_dict[idx]

            # We are in a new block
            if curr_gene != gene_name:

                # 1. Deal with previous gene
                if curr_gene != '' and curr_cds_idx != 0:
                    gene_dict[curr_gene] = [first_start_pos, last_end_pos, curr_cds_start, curr_cds_idx, concat_subseqs]

                # 2. Deal with new gene
                first_start_pos = start_pos
                first_end_pos = end_pos
                curr_gene = gene_name
                curr_cds_idx = 0
                curr_cds_start = -1
                concat_subseqs = dict()
                max_subseq_len = -1
            
            # If we are at the segment with CDS
            if cds_idx != -1:
                curr_cds_idx = max_subseq_len + cds_idx
                curr_cds_start = max_subseq_len
            
            # Add subseqs to concatenated subseq dictionary
            prev_max_subseq_len = max_subseq_len
            for chrom in subseqs:
                if chrom not in concat_subseqs:
                    concat_subseqs[chrom] = ''
                    if prev_max_subseq_len > -1 and first_start_pos != start_pos:
                        concat_subseqs[chrom] = 'N' * prev_max_subseq_len
                concat_subseqs[chrom] = concat_subseqs[chrom] + subseqs[chrom][1].upper()
                this_len = len(concat_subseqs[chrom])
                if this_len > max_subseq_len:
                    max_subseq_len = this_len
            for chrom in concat_subseqs:
                while len(concat_subseqs[chrom]) < max_subseq_len:
                    concat_subseqs[chrom] = concat_subseqs[chrom] + "N"

            last_start_pos = start_pos
            last_end_pos = end_pos
            
        for gene_name in gene_dict:
            out.write("* " + gene_name + "\n")
            out.write("% " + f"{gene_dict[gene_name][0]} {gene_dict[gene_name][1]} {gene_dict[gene_name][2]} {gene_dict[gene_name][3]}\n")
            subseqs = gene_dict[gene_name][4]
            for chrom in subseqs:
                out.write(subseqs[chrom] + "\n")
            out.write("\n")

def process_neg_dict(neg_dict, expanded_coding_regions, ofn):
    with open(ofn, "w") as out:
            
        gene_dict = dict()
        # {gene_a: [leftmost, rightmost, atg_site, {human: actg, chimp: actt, ...}] }

        curr_gene = ''
        for idx in neg_dict:
            start_pos, end_pos, size, gene_name, cds_idx, subseqs = neg_dict[idx]

            # We are in a new block
            if curr_gene != gene_name:

                # 1. Deal with previous gene
                if curr_gene != '':
                    gene_dict[curr_gene] = [first_start_pos, last_end_pos, concat_subseqs]

                # 2. Deal with new gene
                first_start_pos = start_pos
                first_end_pos = end_pos
                curr_gene = gene_name
                concat_subseqs = dict()
                max_subseq_len = -1
            
            # Add subseqs to concatenated subseq dictionary
            prev_max_subseq_len = max_subseq_len
            for chrom in subseqs:
                if chrom not in concat_subseqs:
                    concat_subseqs[chrom] = ''
                    if prev_max_subseq_len > -1 and first_start_pos != start_pos:
                        concat_subseqs[chrom] = 'N' * prev_max_subseq_len
                concat_subseqs[chrom] = concat_subseqs[chrom] + subseqs[chrom][1].upper()
                this_len = len(concat_subseqs[chrom])
                if this_len > max_subseq_len:
                    max_subseq_len = this_len
            for chrom in concat_subseqs:
                while len(concat_subseqs[chrom]) < max_subseq_len:
                    concat_subseqs[chrom] = concat_subseqs[chrom] + "N"

            last_start_pos = start_pos
            last_end_pos = end_pos
            
        for gene_name in gene_dict:
            subseqs = gene_dict[gene_name][2]
            for key in subseqs:
                if "hg38" in key:
                    break
            query = subseqs[key].replace("-", "")
            if "ATG" in query[100:-103]:
                out.write("* " + gene_name + "\n")
                out.write("% " + f"{gene_dict[gene_name][0]} {gene_dict[gene_name][1]} \n")
                for chrom in subseqs:
                    out.write(subseqs[chrom] + "\n")
                out.write("\n")
    print("process_neg_dict done")

def get_corrected_idx(text, query, start_pos, idx, start_codon):
    diff = idx - start_pos
    i = 0
    new_idx = start_pos
    while i < diff or text[new_idx] not in "ACTG":
        if text[new_idx] in "ACTG":
            i += 1
            new_idx += 1
        else:
            new_idx += 1
    while (text[new_idx] not in "ACTG"): new_idx += 1

    return new_idx

def get_start_codon(text, start_pos, idx):
    # idx is idx for whole block
    # see search_idx in get_corrected_idx for local idx
    query = text.replace("-", "")
    query = query.replace("N", "")

    # idx needs to add "-" count

    idx = get_corrected_idx(text, query, start_pos, idx, 'ATG')

    return 'ATG', idx

def conserves_start_codon(text, idx, start_codon):
    query = text[idx:]

    i = 0
    j = 0
    while (j < 3):
        if i == len(query):
            return False
        if query[i] not in 'ACTG':
            i += 1
        else:
            if query[i] != start_codon[j]:
                return False
            i += 1
            j += 1
    return True

def write_labels(text, idx, start_codon, data, counts_dict, key, human):
    text_short = ''
    labels = ''
    omit_species = False
    k = -1
    for j in range(len(text)):
        # Skip "-" and "N" to get local sequence
        if text[j] not in 'ACTG':
            continue
        text_short += text[j]

        # If we have reached the global ATG idx:
        if j == idx:
            if not human and not conserves_start_codon(text, idx, start_codon):
                omit_species = True
                break
            else:
                k = len(labels)
                labels += '1'
                if human: 
                    key += 1
                    counts_dict[key] = 1
                else: counts_dict[key] += 1        # this species conserves ATG
        else:
            labels += '0'

    if not omit_species:
        data.append(text_short[k-100:k+103] + "\n")
        data.append(labels[k-100:k+103] + "\n")

    return counts_dict, key

# Write gene header
def write_header(gene_name, chrom_num, counts_dict, key, start_codon):
    # gene name, chrom #, # conserved, start codon
    header = "* "
    header += gene_name + " "
    header += chrom_num + " "
    header += str(counts_dict[key]) + " "
    header += start_codon + "\n"
    return header

# Process block file into data file
def process_pos_blocks(ifn, ofn, chrom):

    counts_dict = dict()
    key = -1

    with open(ifn, "r") as txt, open(ofn, "w") as out:

        for line in txt:
            if line != "\n": line = line.strip()
            items = line.split(" ")

            # Gene name
            if line.startswith("*"):
                gene_name = items[2].strip('"')         # extract gene name
                chrom_num = chrom.replace("chr", "")    # extract chrom number
                skip_gene = False
                data = []

            # CDS index
            elif line.startswith("%"):
                start_pos = int(items[3])
                idx = int(items[4])
                human = True
            
            # Human sequence
            elif human:
                human = False

                # Get index of ATG site after CDS start
                start_codon, idx = get_start_codon(line, start_pos, idx)
                if start_codon == None: 
                    skip_gene = True
                    
                # Write corresponding labels string
                counts_dict, key = write_labels(line, idx, start_codon, data, counts_dict, key, True)

            # Non-human sequence
            elif line != "\n":
                if skip_gene: continue

                # Write corresponding labels string
                counts_dict, key = write_labels(line, idx, start_codon, data, counts_dict, key, False)

            # End of block
            else:
                if counts_dict[key] < 300:
                    skip_gene = True
                if skip_gene: continue

                header = write_header(gene_name, chrom_num, counts_dict, key, start_codon)
                out.write(header)

                # Write sequences and labels
                for line in data:
                    if len(line) >= 203:
                        out.write(line)
                out.write("\n")

    return counts_dict

def histogram(counts_dict, png_name):
    histo_dict = dict()
    bin_size = 5
    max_size = 500
    for i in range(0, max_size + bin_size, bin_size):
        histo_dict[i] = 0
    for key in counts_dict:
        count = (counts_dict[key]//bin_size) * bin_size
        histo_dict[count] += 1

    # Extract keys and values
    categories = list(range(0, max_size + bin_size, bin_size))  # Bins/categories
    counts = list(histo_dict.values())  # Counts for each bin

    # Create a bar chart
    plt.bar([c + bin_size / 2 for c in categories], counts, width=bin_size, color='pink')

    # Add labels and title
    plt.xlabel('# of start codon matches')
    plt.ylabel('# times that this occurs')
    plt.title(f'Histogram ({len(counts_dict)} genes)')

    # Save the plot
    plt.savefig(png_name)

def process_neg_blocks(ifn, ofn, chrom):

    with open(ifn, "r") as txt, open(ofn, "w") as out:

        for line in txt:
            if line != "\n": line = line.strip()
            items = line.split(" ")

            # Gene name
            if line.startswith("*"):
                gene_name = items[1]                    # extract gene name
                chrom_num = chrom.replace("chr", "")    # extract chrom number
                data = []

            # CDS index
            elif line.startswith("%"):
                human = True
            
            # Human sequence
            elif human:
                human = False

                query = line.replace("-", "").replace("N", "")
                for i in range(100, len(query) - 100):
                    if query[i:i+3] == 'ATG':
                        data.append(query[i-100:i+103])
                        data.append('0'*203)

            # Non-human sequence
            elif line != "\n":

                # Write corresponding labels string
                query = line.replace("-", "").replace("N", "")
                for i in range(100, len(query) - 100):
                    if query[i:i+3] == 'ATG':
                        data.append(query[i-100:i+103]+"\n")
                        data.append('0'*203 + "\n")

            # End of block
            else:
                out.write('* ' + gene_name + "\n")

                # Write sequences and labels
                for line in data:
                    out.write(line)
                out.write("\n")

    print("process_neg_blocks done")

def concatenate_data(pos_data_file, neg_data_file, cat_data, cat_labels):

    with (open(pos_data_file, "r") as pos_txt, open(neg_data_file, "r") as neg_txt,
          open(cat_data, "w") as data, open(cat_labels, "w") as labels):

        pos_lines = pos_txt.readlines()
        neg_lines = neg_txt.readlines()

        # Get positive data sequences
        pos_data = []
        for line in pos_lines:
            if line[0].isalpha():
                pos_data.append(line)

        # Get negative data sequences
        neg_data = []
        for line in neg_lines:
            if line[0].isalpha():
                neg_data.append(line)


        X = len(pos_data)  # Number of lines to sample
        sampled_neg_data = []

        # Sampling algorithm for large files
        for i, line in enumerate(neg_data):
            if i < X:
                sampled_neg_data.append(line)  # Fill reservoir initially
            else:
                j = random.randint(0, i)
                if j < X:
                    sampled_neg_data[j] = line  # Replace an existing sample

        for line in pos_data:
            data.write(line)
        for line in sampled_neg_data:
            data.write(line)

        labels.write("1\n" * len(pos_data))
        labels.write("0\n" * len(sampled_neg_data))

        print(len(pos_data), len(sampled_neg_data))
        

def main(chrom):

    maf_file = "/home/iminkin2/projects3/multiz470way/maf/" + chrom + ".maf"
    gtf_file = "ref/" + chrom + ".mane.gtf"
    block_pos_file = "blocks/pos/" + chrom + ".txt"
    block_neg_file = "blocks/neg/" + chrom + ".txt"
    pos_data_file = "train_data/pos/" + chrom + ".txt"
    neg_data_file = "train_data/neg/" + chrom + ".txt"
    histo_file = "histo/" + chrom + ".png"
    cat_data = "train_data/cat/" + chrom + ".data"
    cat_labels = "train_data/cat/" + chrom + ".labels"

    coding_regions = parse_gtf(gtf_file)
    expanded_coding_regions = get_expanded_regions(coding_regions)
    pos_dict, neg_dict = create_dicts(maf_file, expanded_coding_regions, chrom)
    process_pos_dict(pos_dict, expanded_coding_regions, block_pos_file)
    process_neg_dict(neg_dict, expanded_coding_regions, block_neg_file)

    counts_dict = process_pos_blocks(block_pos_file, pos_data_file, chrom)
    process_neg_blocks(block_neg_file, neg_data_file, chrom)
    histogram(counts_dict, histo_file)

    concatenate_data(pos_data_file, neg_data_file, cat_data, cat_labels)

if __name__ == "__main__":
    main(sys.argv[1])
