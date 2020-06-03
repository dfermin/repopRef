# repopRef
The name is short for _repopulate reference calls_.

This is a tool to "backfill" reference calls in a group of VCF files where they have been dropped.
I had to create this for a project where I was given a bunch of VCF files in which all reference calls were removed.
Each VCF file represented one sample. It was when I tried to merge the VCF files that I found out all the reference calls had been dropped. 

This project requires the [HTSJDK](https://github.com/samtools/htsjdk) which you have to download separately.

When compiled, this JAR will merge multiple VCF files into a single VCF file.
The program assumes each input VCF file represents one sample and only contains non-reference variant calls for that sample.
The program records all the observed calls across all the input VCF files and then writes out a new VCF file to STDOUT.
Where ever a variant call was not observed in sample, it's written out as a reference call. 

So for example consider these 3 snippets from 3 independent VCF files:
```
#CHROM  POS ID  REF ALT INFO  FORMAT sample_1
2 1234  rs1 A T AC=2;AF=1.00,AN=2 GT:AD:DP 1/1:0,1:1
```
```
#CHROM  POS ID  REF ALT INFO  FORMAT sample_2
3 4567  rs2 G C AC=1;AF=1.00,AN=2 GT:AD:DP 0/1:0,1:10
```

```
#CHROM  POS ID  REF ALT INFO  FORMAT sample_3
5 10000  rs3 G C AC=1;AF=1.00,AN=2 GT:AD:DP 0/1:0,1:7
```

When you run the repopRef.jar:
```
java -jar repopRef.jar -i <folder_of_VCFs> | bgzip -c > merged_vcf.gz
```

The output VCF file will look something like this:
```
#CHROM  POS ID  REF ALT INFO  FORMAT sample_1 sample_2  sample_3
2 1234  rs1 A T AC=2;AF=0.33,AN=6 GT:AD:DP 1/1:0,1:1 0/0:-1,-1,1  0/0:-1,-1,1
3 4567  rs2 G C AC=1;AF=0.166,AN=6 GT:AD:DP 0/0:-1,-1:3 0/1:0,1:10  0/0:-1,-1:3
5 10000  rs3 C T AC=1;AF=0.166,AN=6 GT:AD:DP 0/0:-1,-1:2  0/0:-1,-1:2 0/1:0,1:7
```
You will notice that in each row two of the samples have reference calls with `AD` values of `-1,-1` to indicate they are _backfilled_ cases.
For the backfilled reference calls, the `DP` value is the average `DP` value observed among the samples that did have a non-ref call. 


## Warnings
The program loads **all** of the independent VCF files into memory. Depending upon the number and size of your VCF files you may need to break them up by chromosome or region before running them through this program. 
**Memory management is up to you.**


## Licensing
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

