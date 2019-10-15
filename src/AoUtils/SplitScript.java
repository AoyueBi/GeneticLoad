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
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class SplitScript {
    
    public SplitScript(){
        this.splitBwaScript("/Users/Aoyue/Documents/sh_removeBadTaxaFromMergeVCF_notGZ20191015.sh", "sh_removeBadTaxaFromMergeVCF_notGZ_", 10, 4);
        
    }
    
    /**
     * 
     * @param infileS
     * @param nameprefix, the script name you wanna
     * @param numfile, the file number you wanna split to
     * @param numcmd, the number of CDM in each file
     * eg:"/Users/Aoyue/Documents/sh_md5_WheatVMapII_ABgenome_fixmatePosBam.sh", "md5_WheatVMapII_ABgenome_fixmateBam_", 20, 32
     */
    public void splitBwaScript(String infileS,String nameprefix,int numfile, int numcmd) {
        //String infileS = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/001_bwaScript/bwa_20190705needRERUN.sh";
        //String outfileDirS = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/001_bwaScript/splitScript";
        String parentS = new File(infileS).getParent();
        new File(parentS,"splitScript").mkdirs();
        String outfileDirS = new File(parentS,"splitScript").getAbsolutePath();
        String shfileS = new File(parentS,"sh_split.sh").getAbsolutePath();

        try {
            String[] outfileS = new String[numfile];
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter[] bw = new BufferedWriter[numfile];
            for (int i = 0; i < outfileS.length; i++) {
                String num = PStringUtils.getNDigitNumber(3, i + 1);
                outfileS[i] = new File(outfileDirS, nameprefix + num + ".sh").getAbsolutePath();
                bw[i] = IOUtils.getTextWriter(outfileS[i]);
                String temp;
                for (int j = 0; j < numcmd; j++) {
                    if ((temp = br.readLine()) != null) {
                        bw[i].write(temp);
                        bw[i].newLine();
                    }
                }
                bw[i].flush();
                bw[i].close();
            }
            br.close();
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
                bw.write("sh " + fs[i].getName() + " > log_" + fs[i].getName().split(".sh")[0] +".txt &");
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
