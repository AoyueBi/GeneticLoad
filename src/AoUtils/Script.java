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
