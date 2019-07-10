/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vcf;

import java.io.BufferedReader;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.List;
import utils.IOUtils;

/**
 *
 * @author Aoyue
 */
public class CountSites {
    public CountSites(){
        
    }
    
    
    /**
     * Count the snp sites via parallel stream and print it into the inputstream.
     * @param infileDirS 
     */
    public void countSites(String infileDirS){
        //infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_tree/source/";
        File[] fs = new File(infileDirS).listFiles();
        fs=IOUtils.listFilesEndsWith(fs, ".vcf");
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
        System.out.println(new SimpleDateFormat().format(new Date()) + "\tbegin.");
        System.out.println("Chr\tSNP Num");
        fsList.parallelStream().forEach(f -> {
            try{
                String chr = f.getName().split("chr")[1].split(".vcf")[0]; //提取染色体号 001
                int chrint = Integer.parseInt(chr); //将染色体号转化为数字
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp = null;
                int cnt =0;
                while((temp=br.readLine()) != null){
                    if(temp.startsWith("#")) continue;
                    cnt++;
                }
                System.out.println(String.valueOf(chrint) + "\t" + String.valueOf(cnt));
            }
            catch(Exception e){
                e.printStackTrace();
                System.exit(1);
            }
        });
        System.out.println(new SimpleDateFormat().format(new Date()) + "\tend.");
    }
    
    
    public void countVCF(String infileDirS) {
        infileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/originData/";
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        } //注意有隐藏文件，需要进行删除后重新列出文件目录。
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        int sum = 0;
        for (int i = 0; i < fs.length; i++) {
            int cnt = 0;
            try {
                BufferedReader br = IOUtils.getTextGzipReader(fs[i].getAbsolutePath());
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        continue;
                    }
                    cnt++;
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println(String.valueOf(cnt) + "\t" + fs[i].getName());
            sum += cnt;
        }
        System.out.println(String.valueOf(sum));
    }
}
