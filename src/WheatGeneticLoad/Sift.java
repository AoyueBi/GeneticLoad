/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGeneticLoad;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.List;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class Sift {

    public Sift() {
        //this.addGeneBiotype();
        //new FastaBit("/Users/Aoyue/Documents/chr036.fa.gz").getName(0);
        //this.parseGff3();
        //this.mkInputDirs();
        //this.mkParameterConfigChr1_42();
        this.mkJavaCmdchr1_42_SIFTdb();

    }

    public void mkJavaCmdchr1_42_SIFTdb() {
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/001_buildWheatSIFTdb/script/";
        String shfileS = new File("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/001_buildWheatSIFTdb/","sh_mkSIFTdb.sh").getAbsolutePath();
        //nohup java -Xms200g -Xmx500g -jar FastCall.jar parameters_001_FastCall.txt > log_001_fastcall.txt &
        try {

            for (int i = 1; i < 43; i++) {
                String chr = PStringUtils.getNDigitNumber(3, i);
                String outfileS = new File(outfileDirS, "chr" + chr + "_create_wheat_sift4g_db.sh").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write("date\n"
                        + "echo \"chr" + chr + "    :Create wheat sif4G database begins\"" + "\n"
                        + "perl make-SIFT-db-all.pl -config /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/config/abd_iwgscV1_SIFT4G_Create_Genomic_DB_config" + chr + ".txt > log_chr" + chr + ".wheat_sift4G_pipeline.step.txt\n"
                        + "date\n"
                        + "echo \"chr" + chr + "    :Create wheat sif4G database ends\"" + "\n");
                
                bw.flush();
                bw.close();
            }
            

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        
        try {
            File[] fs = new File(outfileDirS).listFiles();
            fs = IOUtils.listFilesEndsWith(fs, ".sh");
            Arrays.sort(fs);
            BufferedWriter bw = IOUtils.getTextWriter(shfileS);
            for (int i = 0; i < fs.length; i++) {
                bw.write("sh " + fs[i].getName() + " > alog_" + fs[i].getName() + " &");
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void mkParameterConfigChr1_42() {
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/001_buildWheatSIFTdb/config/";
        try {
            for (int i = 1; i < 43; i++) {
                String chr = PStringUtils.getNDigitNumber(3, i);
                String outfileS = new File(outfileDirS, "abd_iwgscV1_SIFT4G_Create_Genomic_DB_config" + chr + ".txt").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                StringBuilder sb = new StringBuilder();
                sb.append("\n");
                sb.append("GENETIC_CODE_TABLE=1\n"
                        + "GENETIC_CODE_TABLENAME=Standard\n"
                        + "MITO_GENETIC_CODE_TABLE=11\n"
                        + "MITO_GENETIC_CODE_TABLENAME=Plant Plastids\n").append("\n");
                sb.append("PARENT_DIR=/data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/input").append(chr).append("\n");
                sb.append("ORG=bread_wheat\n"
                        + "ORG_VERSION=abd_iwgscV1\n"
                        + "DBSNP_VCF_FILE=abd_iwgscV1.vcf.gz\n").append("\n");
                sb.append("#Running SIFT 4G\n"
                        + "SIFT4G_PATH=/data1/home/aoyue/biosoftware/SIFT4G/sift4g/bin/sift4g\n"
                        + "PROTEIN_DB=/data1/publicData/protein_db/uniprot90_Jul2018/uniref90.fasta\n"
                        + "COMPUTER=GIS-KATNISS\n"
                        + "\n"
                        + "\n"
                        + "# Sub-directories, don't need to change\n"
                        + "GENE_DOWNLOAD_DEST=gene-annotation-src\n"
                        + "CHR_DOWNLOAD_DEST=chr-src\n"
                        + "LOGFILE=Log.txt\n"
                        + "ZLOGFILE=Log2.txt\n"
                        + "FASTA_DIR=fasta\n"
                        + "SUBST_DIR=subst\n"
                        + "ALIGN_DIR=SIFT_alignments\n"
                        + "SIFT_SCORE_DIR=SIFT_predictions\n"
                        + "SINGLE_REC_BY_CHR_DIR=singleRecords\n"
                        + "SINGLE_REC_WITH_SIFTSCORE_DIR=singleRecords_with_scores\n"
                        + "DBSNP_DIR=dbSNP\n"
                        + "\n"
                        + "# Doesn't need to change\n"
                        + "FASTA_LOG=fasta.log\n"
                        + "INVALID_LOG=invalid.log\n"
                        + "PEPTIDE_LOG=peptide.log\n"
                        + "ENS_PATTERN=ENS\n"
                        + "SINGLE_RECORD_PATTERN=:change:_aa1valid_dbsnp.singleRecord\n"
                        + "");
                bw.write(sb.toString());
                bw.newLine();
                bw.flush();
                bw.close();

            }

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void mkInputDirs() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/001_buildWheatSIFTdb/003_inputFileList/";
        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            File f = new File(infileDirS, "input" + chr);
            f.mkdirs();
            File chrsrc = new File(f, "chr-src");
            chrsrc.mkdirs();
            File geneAnnotationSrc = new File(f, "gene-annotation-src");
            geneAnnotationSrc.mkdirs();
            File dbSNP = new File(f, "dbSNP");
            dbSNP.mkdirs();

            //输出命令的脚本
            try {
                String scriptS = new File(infileDirS, chr + ".sh").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(scriptS);
//                bw.write("cp -f /Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/001_buildWheatSIFTdb/002_splitGTF/chr" + chr +
//                        ".wheat_v1.1_Lulab.addBiotype.gtf.gz " 
//                + geneAnnotationSrc.getAbsolutePath());

                bw.write("cp -f /data1/publicData/wheat/reference/v1.0/byChr/chr" + chr
                        + ".fa.gz /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_d/input" + chr + "/chr-src/ &");
                bw.newLine();
                bw.flush();
                bw.close();

                StringBuilder sb = new StringBuilder("sh ");
                sb.append(scriptS);
                Runtime run = Runtime.getRuntime();
                Process p = run.exec(sb.toString());
                BufferedReader br = new BufferedReader(new InputStreamReader(p.getErrorStream()));
                StringBuilder ssb = new StringBuilder(); //有时候需要将屏幕信息输出到一个文件中，此时建立输出文本路径。
                String line = null;
                while ((line = br.readLine()) != null) {
                    System.out.println(line);
                    ssb.append(line + "\n");
                }
                String result = ssb.toString(); //有必要时输出出来
                p.waitFor();
                System.out.println(sb.toString());
                System.out.println(result);

                //new File(scriptS).delete();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

            //cat /Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_chrList/chrABDgenomeList.txt | awk '{print $1}' | while read a;do echo "cp -f /data1/publicData/wheat/reference/v1.0/byChr/chr${a}.fa.gz /data4/home/aoyue/vmap2/analysis/008_sift/abd_iwgscV1_SIFT4G_db/input${a}/chr-src/ &";done
        }

    }

    /**
     * 将gtf文件拆分成单条染色体
     */
    public void parseGff3() {
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/001_buildWheatSIFTdb/001_gtf/wheat_v1.1_Lulab.addBiotype.gtf.gz";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/008_sift/001_buildWheatSIFTdb/002_splitGTF/";
        String[] outfileS = new String[43];
        //String OutmtFileS = "";
        //"/Users/Aoyue/Documents/Zea_mays.AGPv4.38.gtf";
        for (int i = 0; i < outfileS.length; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            outfileS[i] = new File(outfileDirS, "chr" + chr + ".wheat_v1.1_Lulab.addBiotype.gtf.gz").getAbsolutePath();
            System.out.println(outfileS[i] + "\n");
        }
        try {
            for (int i = 0; i < outfileS.length; i++) {
                String chr = Integer.toString(i);
                BufferedReader br = IOUtils.getTextGzipReader(infileS);
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS[i]);
                String tem = null;
                List<String> l = null;
                while ((tem = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(tem);
                    if (l.get(0).equals(chr)) {
                        StringBuilder sb = new StringBuilder();
                        sb.append(tem);
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
                br.close();
                bw.flush();
                bw.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * modify gtf file, add gene_biotype to the 9th colum
     */
    public void addGeneBiotype() {
        String InputFileS = "/Users/Aoyue/Documents/wheat_v1.1_Lulab.gtf";
        String OutputFileS = "/Users/Aoyue/Documents/wheat_v1.1_Lulab.addBiotype.gtf";
        try {
            BufferedReader br = IOUtils.getTextReader(InputFileS);
            BufferedWriter bw = IOUtils.getTextWriter(OutputFileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                //String plate = null;
                //plate = PStringUtils.fastSplit(temp).get(2);
                //if(plate.equalsIgnoreCase("exon") | plate.equalsIgnoreCase("CDS") | plate.equalsIgnoreCase("stop_codon") | plate.equalsIgnoreCase("start_codon")){
                List<String> l = null;
                l = PStringUtils.fastSplit(temp);
                StringBuilder sb = new StringBuilder();
                sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(2)).append("\t").append(l.get(3)).append("\t").
                        append(l.get(4)).append("\t").append(l.get(5)).append("\t").
                        append(l.get(6)).append("\t").append(l.get(7)).append("\t").append(l.get(8)).append(" ")
                        .append("gene_biotype").append(" ").append("\"protein_coding\";");
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
