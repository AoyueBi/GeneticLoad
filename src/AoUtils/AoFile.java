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
import java.util.List;

/**
 *
 * @author Aoyue
 */
public class AoFile {
    public AoFile(){
        
    }

    /**
     *
     * @param infileS
     */
    public void readheader(String infileS){
        try {
            BufferedReader br = null;
            if (infileS.endsWith(".txt")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".txt.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }else if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
            }
            else if (infileS.endsWith(".xls")) {
                br = IOUtils.getTextReader(infileS);
            }
            String temp = br.readLine();
            List<String> l = PStringUtils.fastSplit(temp);
            int cnt = -1;
            for (int i = 0; i < l.size(); i++) {
                System.out.println(String.valueOf(i) + "\t" + String.valueOf(l.get(i)));
            }

            br.close();

        } catch (Exception e) {
            e.printStackTrace();
        }

    }


    /**
     * 复制某一目录下的所有文件夹，不复制文件
     *
     * @param infileDirS
     * @param outfileDirS
     */
    public void copyFileDirS(String infileDirS, String outfileDirS){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/002_changeChrPos";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/103_snpClassify/003_mergebySub";
        
        
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            System.out.println(fs[i]);
            new File(outfileDirS, fs[i].getName()).mkdirs();
        }
    }
}
