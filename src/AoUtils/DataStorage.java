/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;


import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author Aoyue
 */
public class DataStorage {

    public DataStorage() {
        
//        this.listAllFiles();
//        this.listMD5files();
//        this.getIDseqsubsamples();
//        this.mkGSADB();
//        this.mergeTxt();
        this.countCaseInGroup();

        

    }
    
    /**
     * 
     * 
     */
    public void countCaseInGroup(){
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/006_subfile/LuLab4T_18_subfileList.txt";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/015_addSIFTgroup/chr_D.SNP_anno_addSIFTgroup.txt.gz";
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";
        RowTable<String> t = new RowTable<>(infileS);
        List<String> l = t.getColumn(2);
        System.out.println(l.size() + " list个数");
        Set<String> s = new HashSet<>(l);
        System.out.println(s.size() + " set个数");
        System.out.println(s);
        for(String a : s){
            System.out.println(a + "    " + Collections.frequency(l, a));
        }
    }
    
    /**
     * 目的：将所有sublist文件合并成一个文件
     *
     */
    public void mergeTxt() {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/006_subfile";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/008_GSAdb/subfileList.txt";
        
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/005_MD5";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/008_GSAdb/All_md5.txt";
        
        String infileDirS = "/Users/Aoyue/Documents/IGDB/01_415Lab/007_LabHDD_management/001_HDDstorageList/test";
        String outfileS = "/Users/Aoyue/Documents/IGDB/01_415Lab/007_LabHDD_management/001_HDDstorageList/find.txt";
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        Arrays.sort(fs);
        try {
            String infileS = null;
            BufferedReader br = null;
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            int cnttotal = 0 ; 
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                br = IOUtils.getTextReader(infileS);
                String temp = null; //read header
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    cnttotal++;
                    StringBuilder sb = new StringBuilder();
                    sb.append(temp);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println(String.valueOf(fs[i].getName()) + "\t" + cnt);
            }
            System.out.println(cnttotal + "\t runs in vmap2 project");
            br.close();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        
        /**
         * 
         * LuLab3T_76_subfileList.txt	320
            LuLab3T_77_subfileList.txt	304
            LuLab3T_78_subfileList.txt	304
            LuLab3T_79_subfileList.txt	304
            LuLab3T_80_subfileList.txt	320
            LuLab3T_81_subfileList.txt	320
            LuLab4T_01_subfileList.txt	384
            LuLab4T_02_subfileList.txt	304
            LuLab4T_03_CleanData_15Samples_subfileList.txt	240
            LuLab4T_04_JiaCe_36-Samples_subfileList.txt	288
            LuLab4T_04_subfileList.txt	16
            LuLab4T_16_subfileList.txt	288
            LuLab4T_17_subfileList.txt	320
            LuLab4T_18_subfileList.txt	320
            LuLab4T_21_subfileList.txt	432
            a_supplement_rawData_subfileList.txt	32
            4496	 runs in vmap2 project
         */
    }

    public void mkGSADB() {
//        String inS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/006_subfile/JiaCe_36-Samples_subfileList.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/008_GSAdb/kkkkkkkk.txt";
        
        String inS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/008_GSAdb/subfileList.txt";
//        String inS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/006_subfile/LuLab4T_16_subfileList.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/008_GSAdb/GSAdatabase.txt";

        
        String outtempS = outfileS.replaceFirst(".txt", ".temp.txt");
        try {
            BufferedReader br = IOUtils.getTextReader(inS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            List<String> l = new ArrayList<>();
            String temp = null;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String sm = l.get(0);
                String subsm = l.get(1);
                String experiment = this.getExperimentID(sm);
                String file1 = this.getPairedName(subsm)[0];
                String file2 = this.getPairedName(subsm)[1];
//                String md51 = this.getMD5(file1);
//                String md52 = this.getMD5(file2);
                String md51 = this.getMD5_2(file1);
                String md52 = this.getMD5_2(file2);
                
                bw.write(sm + "\t" + subsm + "\t" + experiment + "\t" + file1 + "\t"+ md51+ "\t" + file2+ "\t" + md52);
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

    public String getMD5_2(String subsm) {
//        String inS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/005_MD5/LuLab4T_04_md5.txt";
        String inS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/008_GSAdb/All_md5.txt";
        String temp = null;
        String md5 = null;
        String out = null;
        List<String> lsubnames = new ArrayList<>();
        HashMap<String,String> hm = new HashMap<>();
        try {
            BufferedReader br = IOUtils.getTextReader(inS);
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplitOnWhitespace(temp);
                md5 = l.get(0);
                String sub = l.get(1);
                if(sub.startsWith("./")){
                    sub = sub.replaceFirst("./", "");
                }
                sub = sub.split("/")[1];
                lsubnames.add(sub);
                hm.put(sub,md5);
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        Collections.sort(lsubnames);
        int index = Collections.binarySearch(lsubnames,subsm);
        if(index > -1){
            md5 = hm.get(lsubnames.get(index));
            out= md5;
        }
        else{
            out = "NA";
            System.out.println("This subsm " + subsm + " has not been found in md5 DataBase");
        }
        
        return out;
    }

    public String getMD5(String subsm) {
        String inS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/005_MD5/LuLab4T_04_md5.txt";
        String temp = null;
        String md5 = null;
        try {
            BufferedReader br = IOUtils.getTextReader(inS);
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplitOnWhitespace(temp);
                String sub = l.get(1);
                if (sub.contains(subsm)) {
                    md5 = l.get(0);
                    return md5;
                }
                else{
                    System.out.println("This subsm" + subsm + "has not been found in md5 DataBase");
                }
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        return md5;

    }

    public String[] getPairedName(String subsm) {
        String[] names = new String[2];
        names[0] = subsm + "_1.clean.fq.gz";
        names[1] = subsm + "_2.clean.fq.gz";
        return names;
    }

    public String getExperimentID(String sm) {
        String inS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/007_HashMapforExperimentIDandSM/ExperimentIDandSM.txt";
        HashMap<String, String> hm = new HashMap();
        RowTable<String> t = new RowTable<>(inS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            String experimentID = t.getCell(i, 0);
            String SM = t.getCell(i, 1);
            hm.put(SM, experimentID);
        }

        return hm.get(sm);
    }

    /**
     * 根据每个SM列出其子样品清单
     *
     */
    public void getIDseqsubsamples() {
//        String infileDirS = "/Volumes/LuLab4T_04/JiaCe_36-Samples"; //必须是包含SAMPLE的目录
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/006_subfile";
//        String infileDirS = "/Volumes/LuLab4T_04";
//        String infileDirS = "/Volumes/LuLab4T_01";
//        String infileDirS = "/Volumes/LuLab4T_02";
//        String infileDirS = "/Volumes/LuLab4T_03/CleanData_15Samples";
//        String infileDirS = "/Volumes/LuLab3T_76";
//        String infileDirS = "/Volumes/LuLab3T_77";
//        String infileDirS = "/Volumes/LuLab3T_78";
//        String infileDirS = "/Volumes/LuLab3T_79";
//        String infileDirS = "/Volumes/LuLab3T_80";
//        String infileDirS = "/Volumes/LuLab3T_81";
//        String infileDirS = "/Volumes/LuLab4T_21";
//        String infileDirS = "/Volumes/LuLab4T_16";
//        String infileDirS = "/Volumes/LuLab4T_17";
//        String infileDirS = "/Volumes/LuLab4T_18";
        String infileDirS = "/Volumes/LuLab4T_18/a_supplement_rawData";
//        String infileDirS = "";
//        String infileDirS = "";
//        String infileDirS = "";
//        String infileDirS = "";
//        String infileDirS = "";
//        String infileDirS = "";
//        String infileDirS = "";
//        String infileDirS = "";
//        String infileDirS = "";
//        String infileDirS = "";
//        String infileDirS = "";
//        

        /**
         * ****************** 建立一个id对应的所有fq文件的数据库 ********************
         */
        try {
            File[] fs = new File(infileDirS).listFiles();
            String infileS = new File(infileDirS).getName();
            String outfileS = new File(outfileDirS, infileS + "_subfileList.txt").getAbsolutePath(); //根据硬盘名字命名
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < fs.length; i++) {
                if (!(fs[i].getName().startsWith("A") || fs[i].getName().startsWith("W") || fs[i].getName().startsWith("B"))) {
                    continue;
                }
                if (fs[i].isDirectory()) {
                    String id = fs[i].getName();
                    File[] fsamples = new File(fs[i].getAbsolutePath()).listFiles();
                    for (int j = 0; j < fsamples.length; j++) {
                        if (fsamples[j].isHidden()) {
                            continue;
                        }
                        if (fsamples[j].isDirectory()) {
                            String sample = fsamples[j].getName();
                            bw.write(id + "\t" + sample + "\n");
                        }
                    }
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
   

    /**
     * 列出以MD5结尾的文件，并将该文件读到一个特定的文件中
     *
     */
    public void listMD5files() {
        double size = Double.MIN_VALUE;
//        String infileDirS = "/Volumes/LuLab4T_04";
//        String infileDirS = "/Volumes/LuLab4T_01";
//          String infileDirS = "/Volumes/LuLab4T_02";
//        String infileDirS = "/Volumes/LuLab4T_03";
//        String infileDirS = "/Volumes/LuLab3T_76";
//        String infileDirS = "/Volumes/LuLab3T_77";
//        String infileDirS = "/Volumes/LuLab3T_78";
//        String infileDirS = "/Volumes/LuLab3T_79";
//        String infileDirS = "/Volumes/LuLab3T_80";
//        String infileDirS = "/Volumes/LuLab3T_81";
//        String infileDirS = "/Volumes/LuLab4T_21";
//        String infileDirS = "/Volumes/LuLab4T_16";
//        String infileDirS = "/Volumes/LuLab4T_17";
        String infileDirS = "/Volumes/LuLab4T_18";
//        String infileDirS = "";
//        String infileDirS = "";
//        String infileDirS = "";
//        String infileDirS = "";
//        String infileDirS = "";
//        String infileDirS = "";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/005_MD5";
        List<File> fsList = new ArrayList();
        File[] fs = null;
        File[] fsparent = IOUtils.listRecursiveFiles(new File(infileDirS)); //列出硬盘下的第一级目录的文件列表，如果是隐藏文件，就跳过不处理
        fs = IOUtils.listFilesEndsWith(fsparent, "md5");

        try {
            String infileS = new File(infileDirS).getName();
            String outfileS = new File(outfileDirS, infileS + "_md5.txt").getAbsolutePath(); //根据硬盘名字命名
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            int cnt = 0;
            int cntmd5 = 0 ;
            for (int i = 0; i < fs.length; i++) {
                if (fs[i].getParentFile().isHidden()) {
                    continue;
                }
                cnt++;
                BufferedReader br = IOUtils.getTextReader(fs[i].getAbsolutePath());
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    bw.write(temp);
                    bw.newLine();
                    cntmd5++;
                }
                br.close();
                System.out.println(fs[i] + " " + cnt + " files has been written to " + outfileS);
            }
            System.out.println(cntmd5/2 + " subsm");
            System.out.println();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    private void listAllFiles() {
        double size = Double.MIN_VALUE;
//        String infileDirS = "/Volumes/LuLab3T_22";
//        String infileDirS = "/Volumes/LuLab3T_23";
//        String infileDirS = "/Volumes/Wheat_Yao"; //截止到硬盘的名字这一级
//        String infileDirS = "/Volumes/Wheat300";
//        String infileDirS = "/Volumes/WheatRef_V3";
//        String infileDirS = "/Volumes/LuLab_wheat";
//        String infileDirS = "/Volumes/lulab45_CS5";
//        String infileDirS = "/Volumes/Lulab3T_46";
//        String infileDirS = "/Volumes/Lulab_47";
//        String infileDirS = "/Volumes/Lulab3T_48";
//        String infileDirS = "/Volumes/LuLab3T_69";
//        String infileDirS = "/Volumes/LuLab3T_70";

//        String infileDirS = "/Volumes/WheatGBS_backup";
//        String infileDirS = "/Volumes/WheatRef";
//        String infileDirS = "/Volumes/LuLab4T_14";
//        String infileDirS = "/Volumes/LuLab4T_23";
//        String infileDirS = "/Volumes/LuLab4T_24";
//        String infileDirS = "/Volumes/LuLab3T_74";
//        String infileDirS = "/Volumes/LuLab3T_72";
//        String infileDirS = "/Volumes/LuLab3T_42";
//        String infileDirS = "/Users/Aoyue/Documents";
//        String infileDirS = "/Volumes/LuLab4T_04";
//        String infileDirS = "/Volumes/LuLab4T_03";
//        String infileDirS = "/Volumes/LuLab3T_76";
//        String infileDirS = "/Volumes/LuLab3T_77";
//        String infileDirS = "/Volumes/LuLab3T_78";
//        String infileDirS = "/Volumes/LuLab3T_79";
//        String infileDirS = "/Volumes/LuLab3T_80";
//        String infileDirS = "/Volumes/LuLab3T_81";
//        String infileDirS = "/Volumes/LuLab4T_21";
//        String infileDirS = "/Volumes/LuLab4T_16";
//        String infileDirS = "/Volumes/LuLab4T_17";
        String infileDirS = "/Volumes/LuLab4T_18";
//        String infileDirS = "";
//        String infileDirS = "";
        String outfileDirS = "/Users/Aoyue/Documents/IGDB/01_415Lab/007_LabHDD_management/001_HDDstorageList";
        List<File> fsList = new ArrayList();
        File[] fs = null;
        File[] fsparent = new File(infileDirS).listFiles(); //列出硬盘下的第一级目录的文件列表，如果是隐藏文件，就跳过不处理
        for (int i = 0; i < fsparent.length; i++) {
            if (fsparent[i].isHidden()) {
                continue;
            }
            if (fsparent[i].getName().startsWith("$") || fsparent[i].getName().startsWith("Sys") || fsparent[i].getName().startsWith("~ ")) {
                continue;
            }
            if (fsparent[i].isFile()) {
                fsList.add(fsparent[i]);
            }
            if (fsparent[i].isDirectory()) {
                fs = IOUtils.listRecursiveFiles(fsparent[i]);
                for (int j = 0; j < fs.length; j++) {
                    if (fs[j].isHidden()) {
                        continue;
                    }
                    fsList.add(fs[j]); //我们最终想要的文件列表
                }
            }
        }

        fs = fsList.toArray(new File[fsList.size()]);
        try {
            String infileS = new File(infileDirS).getName();
            String outfileS = new File(outfileDirS, infileS + ".txt").getAbsolutePath(); //根据硬盘名字命名
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("FileName\tSize(M)");
            bw.newLine();
            for (int i = 0; i < fs.length; i++) {
                if (fs[i].getParentFile().isHidden()) {
                    continue;
                }
                double filesize = fs[i].length() / 1000 / 1000;
                size = size + filesize;
                bw.write(fs[i].getAbsolutePath() + "\t" + String.format("%.0f", filesize));
                bw.newLine();
                System.out.println(fs[i] + " " + String.format("%.0f", filesize));
            }
            System.out.println(String.format("%.0f", size) + "M\tused\t" + String.format("%.0f", (2795000 - size)) + "M\tremained");
            bw.write(String.format("%.0f", size) + "M\tused\t" + String.format("%.0f", (2795000 - size)) + "M\tremained");
            bw.newLine();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
