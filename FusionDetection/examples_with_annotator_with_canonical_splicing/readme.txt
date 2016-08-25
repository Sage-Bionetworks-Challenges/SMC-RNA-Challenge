1.The fusionBedpeAnnotator code in ./fusionBedpeAnnotator:
Compile:
mkdir fusionBedpeAnnotator_build
cd fusionBedpeAnnotator_build
cmake ../fusionBedpeAnnotator/ -DCMAKE_BUILD_TYPE=release 
make
cp ./bin/fusionBedpeAnnotator ../
cd ..

2.Some files needed to run (You should have these already):
(a) GRCh37 fasta
Can be downloaded at ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz.
(b) Annotation file
Run the following commands:
(i)gtfToGenePred -genePredExt -geneNameAsName2 Homo_sapiens.GRCh37.75.gtf Homo_sapiens.GRCh37.75.genePred
(ii) cut -f 1-10,12 Homo_sapiens.GRCh37.75.genePred > tmp.txt
(iii) echo -e "#GRCh37.ensGene.name\tGRCh37.ensGene.chrom\tGRCh37.ensGene.strand\tGRCh37.ensGene.txStart\tGRCh37.ensGene.txEnd\tGRCh37.ensGene.cdsStart\tGRCh37.ensGene.cdsEnd\tGRCh37.ensGene.exonCount\tGRCh37.ensGene.exonStarts\tGRCh37.ensGene.exonEnds\tGRCh37.ensemblToGeneName.value" > annot.enseml.GRCh37.txt
(iv) cat tmp.txt >> annot.enseml.GRCh37.txt

3.The evaluation (assume validated, which is the case for the example res.bedpe):
(a)./fusionBedpeAnnotator -r Homo_sapiens.GRCh37.75.gtf -g annot.enseml.GRCh37.txt -d difile.txt -i truth.bedpe -o truth.annot.bedpe (This only need to run once for one truth file)
(b)./fusionBedpeAnnotator -r Homo_sapiens.GRCh37.75.gtf -g annot.enseml.GRCh37.txt -d difile.txt -i res.bedpe -o res.annot.bedpe
(c)fusionToolEvaluator -t ./truth.annot.bedpe -r res.annot.bedpe -g annot.enseml.GRCh37.txt -s ./rulefile.txt -o evaluate.txt -u

Note: rulefile.examples.txt contains more example rules.
      The Overal, Exon boundary, Canonical splicing, and Exon boundary && Canonical splicing are copied to rulefile.txt already

      You should have fusionToolEvaluator installed already.



