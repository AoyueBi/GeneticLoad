/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import java.io.File;

/**
 *
 * @author Aoyue
 */
public class AoFile {
    public AoFile(){
        
    }
    
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
