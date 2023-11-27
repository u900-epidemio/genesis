#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

# TODO Add additional regexes for new tools in process get_software_versions
regexes = {
    'Pipeline': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'TrimGalore': ['v_trimgalore.txt', r"version (\S+)"],
    'BWA': ['v_bwa.txt', r"Version: (\S+)"],
    'Samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'Bedtools': ['v_bedtools.txt', r"bedtools v(\S+)"],
    'GATK': ['v_gatk.txt', r"Version:(\S+)"],
    'Bcftools': ['v_bcftools.txt', r"bcftools (\S+)"],
    'SnpEff': ['v_snpeff.txt', r"SnpEff version SnpEff (\S+)"],
    'SnpSift': ['v_snpsift.txt', r"SnpSift version (\S+)"],
    'VcfAnno': ['v_vcfanno.txt', r"vcfanno version (\S+)"]
}


results = OrderedDict()
results['Pipeline'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['TrimGalore'] = '<span style="color:#999999;\">N/A</span>'
results['BWA'] = '<span style="color:#999999;\">N/A</span>'
results['Samtools'] = '<span style="color:#999999;\">N/A</span>'
results['Bedtools'] = '<span style="color:#999999;\">N/A</span>'
results['GATK'] = '<span style="color:#999999;\">N/A</span>'
results['Bcftools'] = '<span style="color:#999999;\">N/A</span>'
results['SnpEff'] = '<span style="color:#999999;\">N/A</span>'
results['SnpSift'] = '<span style="color:#999999;\">N/A</span>'
results['VcfAnno'] = '<span style="color:#999999;\">N/A</span>'



# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'Software Versions'
section_href: ''
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd>{}</dd>".format(k,v))
print ("    </dl>")
