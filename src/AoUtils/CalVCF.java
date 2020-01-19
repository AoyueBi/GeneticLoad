/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import format.table.RowTable;
import gnu.trove.list.array.TDoubleArrayList;
import utils.IOUtils;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

/**
 *
 * @author Aoyue
 */
public class CalVCF {

    public CalVCF() {
//        this.calSiteMeanDepth();
//        this.reheader();
//        this.mkPhylipFormat();


    }

    /**
     * 根据barley的fasta格式，写出phylip格式，连续型的
     */
    public void mkPhylipFormat(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/A_subgenome.fasta.txt";
        String outfileS ="/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/chrA_sub_barley_phylip.txt";

        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            List<String> fasta = new ArrayList<>();
            List<String> l = new ArrayList<>();

            while ((temp = br.readLine()) != null) {
                if (temp.startsWith(">")) {
                    continue;
                }
                else {
                    for (int i = 0; i < temp.length(); i++) {
                        fasta.add(temp.substring(i,i+1));
                    }
                }
            }
            br.close();

            //确定要写几行
            int lines = Integer.MIN_VALUE;
            if(fasta.size()%50 == 0){
                lines = fasta.size()/50;
            }
            else{
                lines = fasta.size()/50 + 1;
            }

            //先写第一行
            int k = -1;
            bw.write("Barley");
            int space = 15 - "Barley".length();
            for (int i = 0; i < space; i++) {
                bw.write(" ");
            }
            for (int i = 0; i < 50; i++) {
                k++;
                bw.write(fasta.get(k));
            }
            bw.newLine();

            for (int i = 1; i < lines; i++) {
                //开始写第二行
                for (int j = 0; j < 65; j++) {
                    if(j<15){
                        bw.write(" ");
                    }
                    else{
                        k++;
                        if(k<fasta.size()){
                            bw.write(fasta.get(k));
                        }
                    }
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println(infileS + " is completed at " + outfileS);

        }
        catch (Exception e) {
            e.printStackTrace();
        }

    }


    /**
     * 改变VCF的taxa名字
     *
     */
    public void reheader(String infileS, String outfileS, String reheaderS){
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/chr.Asubgenome.vcf.gz";
//        String outfileS= "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/chr.Asubgenome_reheader.vcf.gz";
//        String reheaderS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";

        RowTable<String> t = new RowTable<>(reheaderS);
        HashMap<String,String> hm = new HashMap<>();
        for (int i = 0; i < t.getRowNumber() ; i++) {
            String taxa = t.getCell(i,0);
            String taxaID = t.getCell(i,1);
            hm.put(taxa,taxaID);
        }

        try {
            BufferedReader br = IOUtils.getTextGzipReader(infileS);
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            String temp = null;
            List<String> taxaList = new ArrayList<>();
            List<String> l = new ArrayList<>();

            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##")) {
                    bw.write(temp);
                    bw.newLine();
                }
                else if(temp.startsWith("#C")){
                    l = PStringUtils.fastSplit(temp);
                    StringBuilder sb = new StringBuilder();
                    for (int i = 0; i < 9; i++) {
                        sb.append(l.get(i)).append("\t");
                    }
                    sb.deleteCharAt(sb.length() - 1);
                    for (int i = 9; i < l.size(); i++) {
                        String taxa = l.get(i);
                        taxa = hm.get(taxa);
                        sb.append("\t").append(taxa);
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                else if(!temp.startsWith("#")) { //
                    bw.write(temp);
                    bw.newLine();

                }//
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(infileS + " is completed at " + outfileS);

        }
        catch (Exception e) {
            e.printStackTrace();
        }

    }


    public TDoubleArrayList calSNPSitesHeter(List<String> vcf){
        TDoubleArrayList out = new TDoubleArrayList();




        return out;
    }

//vcftools --gzvcf hmp321_agpv4_chr1.vcf.gz --weir-fst-pop China_specific.txt --weir-fst-pop Mixed.txt --out China_specificVSMixed_chr001 &
//vcftools --gzvcf hmp321_agpv4_chr1.vcf.gz --weir-fst-pop Non_stiff_stalk.txt --weir-fst-pop China_specific.txt --fst-window-size 10000 --fst-window-step 2000 --out Non_stiff_stalkVSChina_specific_chr001 &
//vcftools --gzvcf hmp321_agpv4_chr10.vcf.gz --keep China_specific.txt --site-pi --out China_specific_hmp321_agpv4_chr10 &
//vcftools --gzvcf hmp321_agpv4_chr10.vcf.gz --keep China_specific.txt --window-pi 10000 --out China_specific_hmp321_agpv4_chr10 &
//vcftools --site-mean-depth --gzvcf /data4/home/aoyue/vmap2/analysis/013_subsetvcf/singleChr/chr036.Dlineage.maf0.005.bi_subset.vcf.gz &
// vcftools --depth --gzvcf /data4/home/aoyue/vmap2/analysis/013_subsetvcf/singleChr/chr036.Dlineage.maf0.005.bi_subset.vcf.gz --out chr036.Dlineage.maf0.005.bi_subset &
    
    public void calSiteMeanDepth() {
        String infileDirS = "";
        String outfileDirS = "";
        List<Integer> lA = new ArrayList<>();
        List<Integer> lD = new ArrayList<>();
        //先进行D的建立
        int j = 5;
        lD.add(j);
        for (int i = 0; i < 6; i++) {
            j = j + 6;
            lD.add(j);
        }
        int k = 6;
        lD.add(k);
        for (int i = 0; i < 6; i++) {
            k = k + 6;
            lD.add(k);
        }
        //再进行A的建立
        int a = 1;
        lA.add(a);
        for (int i = 0; i < 6; i++) {
            a = a + 6;
            lA.add(a);
        }
        int aa = 2;
        lA.add(aa);
        for (int i = 0; i < 6; i++) {
            aa = aa + 6;
            lA.add(aa);
        }

        Collections.sort(lA);
        Collections.sort(lD);
        
        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            int index = Collections.binarySearch(lD, i);
            int index2 = Collections.binarySearch(lA, i);
            if (index < 0) { //说明是属于AB的
                if (index2 > -1) { //说明是属于A的
                    String mPath = new File(infileDirS, "chr" + chr + ".Alineage.vcf").getAbsolutePath();
                    System.out.println("vcftools --site-mean-depth --gzcvf " + mPath );
                } else { //说明是属于B的
                    String mPath = new File(infileDirS, "chr" + chr + ".Blineage.vcf").getAbsolutePath();
                    System.out.println();
                }
            } else if (index > -1) { //说明是属于D的
                String mPath = new File(infileDirS, "chr" + chr + ".Dlineage.vcf").getAbsolutePath();
                System.out.println();
            }
        }
    }
    
    public void bcftools_merge() {
        String abdFileDirS = "/data4/home/aoyue/vmap2/analysis/012_hapscanner/abd/output/VCF/";
        String abFileDirS = "/data4/home/aoyue/vmap2/analysis/012_hapscanner/ab/output/VCF/";
        String dFileDirS = "/data4/home/aoyue/vmap2/analysis/012_hapscanner/d/output/VCF/";
        String mergedFileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/001_rawMergedVCF/";
        /**
         * pseudo-code: 1.建立3个lineage的list,然后进行循环，判断：在A lineage下合并，依次类推。
         */
        List<Integer> lA = new ArrayList<>();
        List<Integer> lD = new ArrayList<>();
        //先进行D的建立
        int j = 5;
        lD.add(j);
        for (int i = 0; i < 6; i++) {
            j = j + 6;
            lD.add(j);
        }
        int k = 6;
        lD.add(k);
        for (int i = 0; i < 6; i++) {
            k = k + 6;
            lD.add(k);
        }
        //再进行A的建立
        int a = 1;
        lA.add(a);
        for (int i = 0; i < 6; i++) {
            a = a + 6;
            lA.add(a);
        }
        int aa = 2;
        lA.add(aa);
        for (int i = 0; i < 6; i++) {
            aa = aa + 6;
            lA.add(aa);
        }

        Collections.sort(lA);
        Collections.sort(lD);
        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            String abdPath = new File(abdFileDirS, "chr" + chr + ".vcf.gz").getAbsolutePath();
            String abPath = new File(abFileDirS, "chr" + chr + ".vcf.gz").getAbsolutePath();
            String dPath = new File(dFileDirS, "chr" + chr + ".vcf.gz").getAbsolutePath();
            int index = Collections.binarySearch(lD, i);
            int index2 = Collections.binarySearch(lA, i);
            if (index < 0) { //说明是属于AB的
                if (index2 > -1) { //说明是属于A的
                    String mPath = new File(mergedFileDirS, "chr" + chr + ".Alineage.vcf").getAbsolutePath();
                    System.out.println("/data1/programs/bcftools-1.8/bcftools merge -m all --force-samples -f PASS,. --threads 2 " + abdPath + " " + abPath + " -o " + mPath + " &");
                } else { //说明是属于B的
                    String mPath = new File(mergedFileDirS, "chr" + chr + ".Blineage.vcf").getAbsolutePath();
                    System.out.println("/data1/programs/bcftools-1.8/bcftools merge -m all --force-samples -f PASS,. --threads 2 " + abdPath + " " + abPath + " -o " + mPath + " &");
                }
            } else if (index > -1) { //说明是属于D的
                String mPath = new File(mergedFileDirS, "chr" + chr + ".Dlineage.vcf").getAbsolutePath();
                System.out.println("/data1/programs/bcftools-1.8/bcftools merge -m all --force-samples -f PASS,. --threads 2 " + abdPath + " " + dPath + " -o " + mPath + " &");
            }
        }
    }

    public void bgzip_ABD() {
        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            System.out.println("bgzip -@ 6 chr" + chr + ".vcf && tabix -p vcf chr" + chr + ".vcf.gz &");
        }
    }

}
