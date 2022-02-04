#########################################################################
# This script takes two vcf files at a time and merges them into one with:
# -> common and unique metadata as well as variants from both the files
# -> label each variant for its source file
# Usage (terminal) "python MergeCVF.py file1.vcf file2.vcf"
#########################################################################

import sys, re

freebayes = sys.argv[1]
freebayes_vcf = open(freebayes, 'r')
freebayes_vcf_read = freebayes_vcf.readlines()

varscan = sys.argv[2]
varscan_vcf = open(varscan, 'r')
varscan_vcf_read = varscan_vcf.readlines()

# This function joins list of lists into a string object
def list2str(lst):
    merge_data = ''
    for elements in lst:
        merge_data = merge_data+'\t'.join(elements)

    return merge_data


# It extract metadata data from input files
def MetadataMerge(lst):
    infolines =  []
    filterlines = []
    formatlines = []
    otherlines = []
    chromline = []
    for metaline in lst:
        meta_line_pattern = re.match('^#', metaline)    # specifier for metadata
        if meta_line_pattern:
            
            if '##INFO=<ID=' in metaline:
                infolines.append(metaline)
            elif '##FILTER=<ID=' in metaline:
                filterlines.append(metaline)
            elif '##FORMAT=<ID=' in metaline:
                formatlines.append(metaline)
            elif '#CHROM' in metaline:
                chromline.append(metaline)
            else:
                otherlines.append(metaline)

    return otherlines, infolines, filterlines, formatlines, chromline


# Extract meta data for both the vcf files
freebayes_meta = MetadataMerge(freebayes_vcf_read)
varscan_meta = MetadataMerge(varscan_vcf_read)


# Create a common metadata (provided common values as well as unique values
# for each file are represented once)from both the input file
info_merge = freebayes_meta[1]+varscan_meta[1]
info_merge_uni = list(set(info_merge))

filter_merge = freebayes_meta[2]+varscan_meta[2]
filter_merge_uni = list(set(filter_merge))

format_merge = freebayes_meta[3]+varscan_meta[3]
format_merge_uni = list(set(format_merge))

other_str = ''.join(freebayes_meta[0])
info_str = ''.join(info_merge_uni)
filter_str = ''.join(filter_merge_uni)
format_str = ''.join(format_merge_uni)
chrom_str = freebayes_meta[4][0]


# Extract all chromosomal positions from both the input files and store as a list
def UniuqePos(lst):
    chr_pos = []
    for line in lst:    
        var_lines_pattern = re.match('^chr', line)  # specifier for variant data lines        
        if var_lines_pattern:          
            data_line_split = line.split('\t')
            chr_str = data_line_split[0][0:3]
            chr_num = data_line_split[0][3:]
            data_line_pos = data_line_split[1]
            data_chr_pos = chr_str, int(chr_num), int(data_line_pos)
            chr_pos.append(data_chr_pos)

    return chr_pos


# Unique chromosaml position from both the files
freebayes_pos = UniuqePos(freebayes_vcf_read)
varscan_pos = UniuqePos(varscan_vcf_read)
all_pos = freebayes_pos + varscan_pos
all_pos_uni = [list(x) for x in set(tuple(x) for x in all_pos)]
all_pos_uni_sorted = sorted(all_pos_uni, key=lambda x:(x[1],x[2]))


# It extracts variant data line from input file
def vcfMerge(lst2, vcf):
    data_line = []
    for line in vcf:
        var_lines_pattern = re.match('^chr', line)  # specifier
        if var_lines_pattern:
            data_line_split = line.split('\t')
            lst2str = lst2[0]+str(lst2[1]), str(lst2[2])
            if lst2str[0] == data_line_split[0] and lst2str[1] == data_line_split[1]:                
                data_line = line
                
    return data_line


# Look for common and unique variants from both the files and store
# those variants as list. Each variants are labeled according to their source
vcf_merge_lst = []
for entity in all_pos_uni_sorted:
    vcf1match = vcfMerge(entity,freebayes_vcf_read)
    vcf2match = vcfMerge(entity,varscan_vcf_read)
    if vcf1match != [] and vcf2match != []:
        vcf1match_split = vcf1match.split('\t')
        vcf2match_split = vcf2match.split('\t')
        if vcf1match_split[4]==vcf2match_split[4]:
            vcf1match_split[7] = vcf1match_split[7]+';calledBy=Freebayes+VarScan'
            vcf_merge_lst.append(vcf1match_split) 
        else:
            vcf1match_split[7] = vcf1match_split[7]+';calledBy=Freebayes'
            vcf2match_split[7] = vcf2match_split[7]+';calledBy=VarScan'
            vcf_merge_lst.append(vcf1match_split)
            vcf_merge_lst.append(vcf2match_split)

    elif vcf1match != [] and vcf2match == []:
        vcf1match_split = vcf1match.split('\t')
        vcf1match_split[7] = vcf1match_split[7]+';calledBy=Freebayes'
        vcf_merge_lst.append(vcf1match_split)

    elif vcf1match == [] and vcf2match != []:
        vcf2match_split = vcf2match.split('\t')
        vcf2match_split[7] = vcf2match_split[7]+';calledBy=VarScan'
        vcf_merge_lst.append(vcf2match_split)
    

# List to string conversion of each variants
merge_data_var = list2str(vcf_merge_lst)


# Join metadata and variant data into a string
all_merge_str = other_str + info_str + filter_str + format_str + chrom_str + merge_data_var


# Write output in vcf format 
file = open("MergeOutput.vcf", "w")
file.write(all_merge_str)
file.close()
