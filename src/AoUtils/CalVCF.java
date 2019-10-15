/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class CalVCF {

    public CalVCF() {
        this.siteMeanDepth();

    }

//vcftools --gzvcf hmp321_agpv4_chr1.vcf.gz --weir-fst-pop China_specific.txt --weir-fst-pop Mixed.txt --out China_specificVSMixed_chr001 &
//vcftools --gzvcf hmp321_agpv4_chr1.vcf.gz --weir-fst-pop Non_stiff_stalk.txt --weir-fst-pop China_specific.txt --fst-window-size 10000 --fst-window-step 2000 --out Non_stiff_stalkVSChina_specific_chr001 &
//vcftools --gzvcf hmp321_agpv4_chr10.vcf.gz --keep China_specific.txt --site-pi --out China_specific_hmp321_agpv4_chr10 &
//vcftools --gzvcf hmp321_agpv4_chr10.vcf.gz --keep China_specific.txt --window-pi 10000 --out China_specific_hmp321_agpv4_chr10 &
//vcftools --site-mean-depth --gzvcf /data4/home/aoyue/vmap2/analysis/013_subsetvcf/singleChr/chr036.Dlineage.maf0.005.bi_subset.vcf.gz & 
    //vcftools --depth --gzvcf /data4/home/aoyue/vmap2/analysis/013_subsetvcf/singleChr/chr036.Dlineage.maf0.005.bi_subset.vcf.gz --out chr036.Dlineage.maf0.005.bi_subset & 
    
    public void siteMeanDepth() {
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
