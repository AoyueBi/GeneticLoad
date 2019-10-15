/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import utils.IOUtils;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.io.File;
import java.util.*;

/**
 *
 * @author Aoyue
 */
public class Script {

    public Script() {
        //this.removeBadTaxafromVCF();

    }

    /**
     * find -name "*.vcf" | cut -c3- ; 本地获取运行脚本，输出在netbeans的output界面。
     */
    public void removeBadTaxafromVCF() {
        try {
            String infileS = "/Users/Aoyue/Documents/vcflist.txt";
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                String chr = temp.substring(3, 6);
//vcftools --vcf /data2/aoyue/fastcall_ABgenome/003_rawVCF_removeBadTaxa/chr001.ABgenome.10000lines.vcf --remove-indv B0043  --remove-indv B0205 --remove-indv B0092 --remove-indv B0144 --remove-indv B0028 --remove-indv B0191 --remove-indv B0178 --remove-indv B0046 --recode --recode-INFO-all --out chr001.ABgenome.10000lines.removeBadTaxa
                System.out.println("vcftools --gzvcf /data2/aoyue/fastcall_ABgenome/002_bivcf/" + temp
                        + " --remove /data2/aoyue/fastcall_ABgenome/003_rawVCF_removeBadTaxa/removeList_ABgenome.txt --recode --recode-INFO-all --stdout | bgzip -c -@ 4 > "
                        + "/data2/aoyue/fastcall_ABgenome/004_bivcf_removeBadTaxa/" + temp.replaceFirst(".vcf.gz", ".removeBadTaxa.vcf.gz")
                        + " &");
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 1 2 3 4 5 6 后缀分别是Alineage Blineage Dlineage chr024.Dlineage.vcf
     */
    public void bgzip_lineage() {
        int[] arra = {1,2,7,8,13,14,19,20,25,26,31,32,37,38};
        int[] arrb = {3,4,9,10,15,16,21,22,27,28,33,34,39,40};
        int[] arrd = {5,6,11,12,17,18,23,24,29,30,35,36,41,42};
        HashMap<Integer,String> hml = new HashMap<>();
        Arrays.sort(arra);
        Arrays.sort(arrb);
        Arrays.sort(arrd);
        for(int i =0; i<arra.length;i++){
            hml.put(arra[i], "A");
            hml.put(arrb[i], "B");
            hml.put(arrd[i], "D");
        }
        for(int i=0; i < 42;i++){
            int j = i+1;
            String chr = PStringUtils.getNDigitNumber(3, j);
            String lineage = hml.get(j);
            if(lineage.equals("D")){
                System.out.println("vcftools --gzvcf /data4/home/aoyue/vmap2/genotype/mergedVCF/001_rawMergedVCF/chr" + chr + "." + lineage + "lineage.vcf.gz --remove /data4/home/aoyue/vmap2/analysis/017_removeBadTaxa/000_badTaxaList/BadTaxa_Hexa_Diploid_S8.txt --recode --recode-INFO-all --stdout > /data4/home/aoyue/vmap2/genotype/mergedVCF/004_rawMergedVCF_removeBadTaxa/chr" + chr + ".lineage.vcf" );
                
            }
            else{
                System.out.println("vcftools --gzvcf /data4/home/aoyue/vmap2/genotype/mergedVCF/001_rawMergedVCF/chr" + chr + "." + lineage + "lineage.vcf.gz --remove /data4/home/aoyue/vmap2/analysis/017_removeBadTaxa/000_badTaxaList/BadTaxa_Hexa_Tetra_S8.txt --recode --recode-INFO-all --stdout > /data4/home/aoyue/vmap2/genotype/mergedVCF/004_rawMergedVCF_removeBadTaxa/chr" + chr + ".lineage.vcf" );
                
            }
            
            //System.out.println("bgzip -@ 10 chr" + chr + "." + lineage + "lineage.maf0.005.bi.vcf");
            //System.out.println("chr"+chr+"."+lineage+"lineage.maf0.005.bi.vcf");
            // chr002.Alineage.maf0.005.bi_SIFTannotations.xls
            //System.out.println("cp /data4/home/aoyue/vmap2/analysis/008_sift/result/001_annotatorResult/output"+chr+ "/chr"+ chr+"."+lineage+"lineage.maf0.005.bi_SIFTannotations.xls /data4/home/aoyue/vmap2/analysis/008_sift/result/001_annotatorResult/output/ &");
            //System.out.println("mv /data4/home/aoyue/vmap2/analysis/015_annoDB/001_step1/chr"+chr+"."+lineage+"lineage.maf0.005.bi.AnnoDB.txt.gz /data4/home/aoyue/vmap2/analysis/015_annoDB/001_step1/chr" +chr+"."+ "lineage.maf0.005.bi.AnnoDB.txt.gz");

        }
    }

    public void bgzip_AB() {
        List<Integer> l = new ArrayList<>();
        int j = 5;
        l.add(j);
        for (int i = 0; i < 6; i++) {
            j = j + 6;
            l.add(j);
        }

        int k = 6;
        l.add(k);
        for (int i = 0; i < 6; i++) {
            k = k + 6;
            l.add(k);
        }
        Collections.sort(l);

        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            int index = Collections.binarySearch(l, i);
            if (index < 0) {
                System.out.println("bgzip -@ 6 chr" + chr + ".vcf && tabix -p vcf chr" + chr + ".vcf.gz &");
            }
        }
    }

    /**
     * 压缩文件bgzip并建立索引 chr027.Dgenome.vcf
     */
    public void bgzip_D() {
        List<Integer> l = new ArrayList<>();
        int j = 5;
        l.add(j);
        for (int i = 0; i < 6; i++) {
            j = j + 6;
            l.add(j);
        }

        int k = 6;
        l.add(k);
        for (int i = 0; i < 6; i++) {
            k = k + 6;
            l.add(k);
        }
        Collections.sort(l);

        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            int index = Collections.binarySearch(l, i);
            if (index > -1) {
                //System.out.println("bgzip -@ 6 chr" + chr + ".vcf && tabix -p vcf chr" + chr + ".vcf.gz &");
                System.out.println("bgzip -@ 10 chr" + chr + ".Dgenome.vcf");
            }
        }
    }

    /**
     * chr035.ABDgenome.vcf
     */
    public void bgzip_ABD() {
        for (int i = 1; i < 43; i++) {
            String chr = PStringUtils.getNDigitNumber(3, i);
            //System.out.println("bgzip -@ 10 chr" + chr + ".vcf && tabix -p vcf chr" + chr + ".vcf.gz &");
            System.out.println("bgzip -@ 10 chr" + chr + ".ABDgenome.vcf");
        }
    }

    /**
     * find | cut -f2 -d"/" bgzip -c -@ 10 chr005.vcf > chr005.vcf.gz &
     * 如果不写路径的话，会直接压缩覆盖原来的文件;如果写路径，则会重新生成一个文件，原来未压缩的文件依旧存在。前提是bgzip 不加 -c参数
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void bgzip_deprecated(String infileDirS, String outfileDirS, String threads) {
        /**
         * ** need to modify ***
         */
        //===========================
        String cmd = null;
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf");
        Arrays.sort(fs);
        for (int i = 0; i < fs.length; i++) {
            StringBuilder sb = new StringBuilder();
            sb.append("bgzip -c -@ " + threads + " " + new File(infileDirS, fs[i].getName()).getAbsolutePath() + " > " + new File(outfileDirS, fs[i].getName().replaceFirst(".vcf", ".vcf.gz")).getAbsolutePath() + " &");
            cmd = sb.toString();
            System.out.println(cmd);
        }
    }

}
