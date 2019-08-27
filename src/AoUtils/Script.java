/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import utils.IOUtils;

/**
 *
 * @author Aoyue
 */
public class Script {

    public Script() {
        this.removeBadTaxafromVCF();

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
    
    
    
    public void script_local(String infileDirS, String outfileDirS) {
        /**
         * ** need to modify ***
         */
        int threads = 10;
        //String dbfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_chrList/chrABDgenomeList.txt";
        
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
    
    /**
     * find | cut -f2 -d"/"
     * 如果不写路径的话，会直接压缩覆盖原来的文件;如果写路径，则会重新生成一个文件，原来未压缩的文件依旧存在。q前提是bgzip 不加 -c参数 
     * @param infileDirS
     * @param outfileDirS 
     */
    public void bgzip_noscript(String infileDirS, String outfileDirS) {
        /**
         * ** need to modify ***
         */
        int threads = 10;
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

    /**
     * bgzip -c -@ 10 chr005.vcf > chr005.vcf.gz &
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void bgzip(String scriptS, String infileDirS, String outfileDirS) {
        /**
         * ** need to modify ***
         */
        int threads = 10;
        //===========================
        String cmd = null;
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf");
        Arrays.sort(fs);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(scriptS);

            for (int i = 0; i < fs.length; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append("bgzip -c -@ " + threads + " " + new File(infileDirS, fs[i].getName()).getAbsolutePath() + " > " + new File(outfileDirS, fs[i].getName().replaceFirst(".vcf", ".vcf.gz")).getAbsolutePath() + " &");
                cmd = sb.toString();
                bw.write(cmd);
                bw.newLine();
            }

            System.out.println(cmd);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

}
