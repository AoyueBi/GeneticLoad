/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGeneticLoad;

import format.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class Treetest {

    public Treetest() {
//        this.colbyType();
       //this.colbyContinent();
//        this.binarybyType();
//        this.labels();
        //this.TREE_COLORS();
        //this.colRanges();
        //this.colRangesbyDico();
        //this.colRangebyPloidy();
        //this.prune_removeNA();
        //this.colRangebyPloidy_Dsubgenome();
        //this.prune_removeNA_Dsubgenome();
        //this.colRangebyHexaDiGroup_Dsubgenome();
        //this.colRangebyHexaDiGroup_ABsubgenome();

    }
    
    /**
     * 进行A_B subgenome 的分组
     */
    public void colRangebyHexaDiGroup_ABsubgenome() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/004_A_Bsubgenome_maf0.01/source/002_subspecies/";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/004_A_Bsubgenome_maf0.01/003_labels/001_colRangebyploidy.A_Bsubgenome.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/004_A_Bsubgenome_maf0.01/003_labels/001_colBranchbyploidy.A_Bsubgenome.txt";
        File f = new File(infileDirS);
        File[] fs = IOUtils.listRecursiveFiles(f);
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
            //System.out.println(fs[i].getName().split(".txt")[0]);
        }
        fs = IOUtils.listRecursiveFiles(f);
        Arrays.sort(fs);

        String[] groups = {"Ae.tauschii", "Wild emmer","Domesticated emmer","Free-threshing emmer","Not free-threshing emmer","Cultivar","Landrace","Others"};
        String[] col = { "#87cef9","#ffcf66","#cc8b00","#664601","#b17902","#960505","#fc6e6e","#de0707"};
        HashMap<String, String> hm = new HashMap<>();
        for (int i = 0; i < groups.length; i++) {
            hm.put(groups[i], col[i]);
        }
        try {
            BufferedWriter[] bw = new BufferedWriter[2];
            bw[0] = IOUtils.getTextWriter(outfileS);
            bw[1] = IOUtils.getTextWriter(outfileS1);
            //先写表头
            bw[0].write("TREE_COLORS\nSEPARATOR TAB\nDATA\n");
            bw[1].write("TREE_COLORS\nSEPARATOR TAB\nDATA\n");
            //再写内部的分组，注意文件名字必须和HashMap里的分组名保持一致
            for (int i = 0; i < fs.length; i++) {
                String infileS = fs[i].getAbsolutePath();
                String group = fs[i].getName().split(".txt")[0];
                //if(group.equals("NA")) continue;
                BufferedReader br = IOUtils.getTextReader(infileS);
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    bw[0].write(temp + "\trange\t" + hm.get(group) + "\t" + group);
                    bw[0].newLine();
                    
                    bw[1].write(temp + "\tbranch\t" + hm.get(group) + "\t" + "normal");
                    bw[1].newLine();
                }
                br.close();

            }
            bw[0].flush();
            bw[1].flush();
            bw[0].close();
            bw[1].close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * 进行D subgenome 的分组
     */
    public void colRangebyHexaDiGroup_Dsubgenome() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/003_Dsubgenome_maf0.001/source/003_group/";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/003_Dsubgenome_maf0.001/003_labels/001_colRangebyploidy.Dsubgenome.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/003_Dsubgenome_maf0.001/003_labels/001_colBranchbyploidy.Dsubgenome.txt";
        File f = new File(infileDirS);
        File[] fs = IOUtils.listRecursiveFiles(f);
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
            //System.out.println(fs[i].getName().split(".txt")[0]);
        }
        fs = IOUtils.listRecursiveFiles(f);
        Arrays.sort(fs);

        String[] groups = {"Ae.tauschii", "Wild emmer","Domesticated emmer","Free-threshing emmer","Not free-threshing emmer","Cultivar","Landrace","Others"};
        String[] col = { "#87cef9","#ffcf66","#cc8b00","#664601","#b17902","#960505","#fc6e6e","#de0707"};
        HashMap<String, String> hm = new HashMap<>();
        for (int i = 0; i < groups.length; i++) {
            hm.put(groups[i], col[i]);
        }
        try {
            BufferedWriter[] bw = new BufferedWriter[2];
            bw[0] = IOUtils.getTextWriter(outfileS);
            bw[1] = IOUtils.getTextWriter(outfileS1);
            //先写表头
            bw[0].write("TREE_COLORS\nSEPARATOR TAB\nDATA\n");
            bw[1].write("TREE_COLORS\nSEPARATOR TAB\nDATA\n");
            //再写内部的分组，注意文件名字必须和HashMap里的分组名保持一致
            for (int i = 0; i < fs.length; i++) {
                String infileS = fs[i].getAbsolutePath();
                String group = fs[i].getName().split(".txt")[0];
                //if(group.equals("NA")) continue;
                BufferedReader br = IOUtils.getTextReader(infileS);
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    bw[0].write(temp + "\trange\t" + hm.get(group) + "\t" + group);
                    bw[0].newLine();
                    
                    bw[1].write(temp + "\tbranch\t" + hm.get(group) + "\t" + "normal");
                    bw[1].newLine();
                }
                br.close();

            }
            bw[0].flush();
            bw[1].flush();
            bw[0].close();
            bw[1].close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }
    
    
    /**
     * 
     * @param infileS1 小麦参考基因组文件
     * @param infileS2 输入vcf文件
     * @param outfileS 输出fa文件
     */
    public void getvcf2fasta(String infileS1,String infileS2,String outfileS){
        try{
            String temp = null;
            String temp2 = null;
            int i,j,m;          
            BufferedReader brfasta = IOUtils.getTextReader(infileS1);
            //BufferedReader brsnpAndindel = IOUtils.getTextReader(infileS2);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String fastafirst = brfasta.readLine();
            int startcodon = Integer.valueOf(fastafirst.split(":")[1].split("-")[0]); //1:35-150
            StringBuilder fasta = new StringBuilder();
            Set<Integer> pos = new HashSet();
            HashMap<Integer, String> hashMapfasta = new HashMap<Integer, String>();
            while((temp = brfasta.readLine()) != null){
                fasta.append(temp);         
            }
            //List<String> fastalist = new ArrayList<>();
            for(i=0;i<fasta.length();i++){ //fasta文件的内容
                int key = startcodon +i ; //第key个碱基的位置
                String value = fasta.substring(i, i+1); //第key个碱基
                pos.add(key);
                hashMapfasta.put(key, value);
            }
            
            //String[] firstline2 = brsnpAndindel.readLine().split("\t");
            Set<String> namepos = new HashSet();
//            for(j=5;j<firstline2.length;j++){
//                name.add(firstline2[j]);              
//            }
            RowTable<String> t = new RowTable<>(infileS2);
            //String str[][] = new String[3][4];
            //int aa = t.getColumnNumber();
            //System.out.println(aa);
            List<String> aa = t.getColumn(1); //得到vcf文件的pos列表
            for(int bb = 0;bb<aa.size();bb++){
                System.out.println(aa.get(bb));
                namepos.add(aa.get(bb)); //namepos是vcf所有位点的集合
            }          
            for(m=5; m<t.getColumnNumber();m++){  //进行每个taxa的循环，数列读取
                bw.write("\n");
                List<String> l = t.getColumn(m); //得到该列对应的taxa下的所有位点信息
                int qq = 0;
                System.out.println(l);
                bw.write(">" + t.getColumnName(m) + "\n");
                for(int p = startcodon;p < fasta.length() + startcodon;p++){
                    if(namepos.add(String.valueOf(p))){ //如果能加入该set,说明该位点不是变异位点。
                        namepos.remove(String.valueOf(p));
                        bw.write(hashMapfasta.get(p)); //得到该位点对应的碱基信息
                    }
                    else{
                        //int qq = 0;
                        String aaa = t.getCellAsString(qq, m); //第m列的第0行，接下来会有第1行第2行的判断                      
                        if(aaa.equals("0")){
                            System.out.println(t.getCellAsString(qq, 3));
                            bw.write(t.getCellAsString(qq, 3));
                        }
                        else{
                            System.out.println(t.getCellAsString(qq, 4));
                            bw.write(t.getCellAsString(qq, 4));
                        }
                        qq++;
                    }
                }              
            }
            bw.flush();
            bw.close();
            
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    /**
     * 根据倍性不同进行着色区分
     */
    //
    public void prune_removeNA_Dsubgenome() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/002_ABD_D/source/002_subspecies";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/002_ABD_D/002_labels/002_excludeNAlist.txt";
        File f = new File(infileDirS);
        File[] fs = IOUtils.listRecursiveFiles(f);
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
            //System.out.println(fs[i].getName().split(".txt")[0]);
        }
        fs = IOUtils.listRecursiveFiles(f);
        Arrays.sort(fs);

        String[] groups = {"Ae.tauschii", "Wild emmer","Domesticated emmer","Durum","Carthlicum","Ispahanicum",
            "Karamyschevii","Polonicum","Turanicum","Turgidum","Cultivar","Landrace","Breeding_Research Material"};
        String[] col = {"#87cef9", "#ffcf66", "#cc8c00", "#664600", "#996900", "#996900", "#996900", "#996900", "#996900", "#996900", "#e10505", "#fc6e6e", "#af0404"};
        HashMap<String, String> hm = new HashMap<>();
        for (int i = 0; i < groups.length; i++) {
            hm.put(groups[i], col[i]);
        }
        try {
            BufferedWriter[] bw = new BufferedWriter[2];
            bw[0] = IOUtils.getTextWriter(outfileS);
            //再写内部的分组，注意文件名字必须和HashMap里的分组名保持一致
            for (int i = 0; i < fs.length; i++) {
                String infileS = fs[i].getAbsolutePath();
                String group = fs[i].getName().split(".txt")[0];
                if(group.equals("NA")) continue;
                BufferedReader br = IOUtils.getTextReader(infileS);
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    bw[0].write(temp );
                    bw[0].newLine();
                    
                }
                br.close();

            }
            bw[0].flush();
            bw[0].close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }
    
    /**
     * 根据分组文件，进行合并，得到一个完整的分组信息文件
     */
    public void colRangebyPloidy_Dsubgenome() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/002_ABD_D/source/002_subspecies";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/002_ABD_D/002_labels/001_colRangebyploidy.Dsubgenome.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/002_ABD_D/002_labels/001_colBranchbyploidy.Dsubgenome.txt";
        File f = new File(infileDirS);
        File[] fs = IOUtils.listRecursiveFiles(f);
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
            //System.out.println(fs[i].getName().split(".txt")[0]);
        }
        fs = IOUtils.listRecursiveFiles(f);
        Arrays.sort(fs);

        String[] groups = {"Ae.tauschii", "Wild emmer","Domesticated emmer","Durum","Carthlicum","Ispahanicum",
            "Karamyschevii","Polonicum","Turanicum","Turgidum","Cultivar","Landrace","Breeding_Research Material"};
        String[] col = {"#87cef9", "#ffcf66", "#cc8c00", "#664600", "#996900", "#996900", "#996900", "#996900", "#996900", "#996900", "#e10505", "#fc6e6e", "#af0404"};
        HashMap<String, String> hm = new HashMap<>();
        for (int i = 0; i < groups.length; i++) {
            hm.put(groups[i], col[i]);
        }
        try {
            BufferedWriter[] bw = new BufferedWriter[2];
            bw[0] = IOUtils.getTextWriter(outfileS);
            bw[1] = IOUtils.getTextWriter(outfileS1);
            //先写表头
            bw[0].write("TREE_COLORS\nSEPARATOR TAB\nDATA\n");
            bw[1].write("TREE_COLORS\nSEPARATOR TAB\nDATA\n");
            //再写内部的分组，注意文件名字必须和HashMap里的分组名保持一致
            for (int i = 0; i < fs.length; i++) {
                String infileS = fs[i].getAbsolutePath();
                String group = fs[i].getName().split(".txt")[0];
                if(group.equals("NA")) continue;
                BufferedReader br = IOUtils.getTextReader(infileS);
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    bw[0].write(temp + "\trange\t" + hm.get(group) + "\t" + group);
                    bw[0].newLine();
                    
                    bw[1].write(temp + "\tbranch\t" + hm.get(group) + "\t" + "normal");
                    bw[1].newLine();
                }
                br.close();

            }
            bw[0].flush();
            bw[1].flush();
            bw[0].close();
            bw[1].close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }
    
    

    /**
     * 根据倍性不同进行着色区分
     */
    //
    public void prune_removeNA() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/001_ABsubgenome/source/002_subspecies/";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/001_ABsubgenome/003_labels/002_excludeNAList.txt";
        File f = new File(infileDirS);
        File[] fs = IOUtils.listRecursiveFiles(f);
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
            //System.out.println(fs[i].getName().split(".txt")[0]);
        }
        fs = IOUtils.listRecursiveFiles(f);
        Arrays.sort(fs);

        String[] groups = {"Ae.tauschii", "Wild emmer","Domesticated emmer","Durum","Carthlicum","Ispahanicum",
            "Karamyschevii","Polonicum","Turanicum","Turgidum","Cultivar","Landrace","Breeding_Research Material"};
        String[] col = {"#87cef9", "#ffcf66", "#cc8c00", "#664600", "#996900", "#996900", "#996900", "#996900", "#996900", "#996900", "#e10505", "#fc6e6e", "#af0404"};
        HashMap<String, String> hm = new HashMap<>();
        for (int i = 0; i < groups.length; i++) {
            hm.put(groups[i], col[i]);
        }
        try {
            BufferedWriter[] bw = new BufferedWriter[2];
            bw[0] = IOUtils.getTextWriter(outfileS);
            //再写内部的分组，注意文件名字必须和HashMap里的分组名保持一致
            for (int i = 0; i < fs.length; i++) {
                String infileS = fs[i].getAbsolutePath();
                String group = fs[i].getName().split(".txt")[0];
                if(group.equals("NA")) continue;
                BufferedReader br = IOUtils.getTextReader(infileS);
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    bw[0].write(temp );
                    bw[0].newLine();
                    
                }
                br.close();

            }
            bw[0].flush();
            bw[0].close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     * 根据分组文件，进行合并，得到一个完整的分组信息文件
     */
    public void colRangebyPloidy() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/001_ABsubgenome/source/002_subspecies";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/001_ABsubgenome/003_labels/001_colRangebyploidy.ABsubgenome.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/001_ABsubgenome/003_labels/001_colBranchbyploidy.ABsubgenome.txt";
        File f = new File(infileDirS);
        File[] fs = IOUtils.listRecursiveFiles(f);
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
            //System.out.println(fs[i].getName().split(".txt")[0]);
        }
        fs = IOUtils.listRecursiveFiles(f);
        Arrays.sort(fs);

        String[] groups = {"Ae.tauschii", "Wild emmer","Domesticated emmer","Durum","Carthlicum","Ispahanicum",
            "Karamyschevii","Polonicum","Turanicum","Turgidum","Cultivar","Landrace","Breeding_Research Material"};
        String[] col = {"#87cef9", "#ffcf66", "#cc8c00", "#664600", "#996900", "#996900", "#996900", "#996900", "#996900", "#996900", "#e10505", "#fc6e6e", "#af0404"};
        HashMap<String, String> hm = new HashMap<>();
        for (int i = 0; i < groups.length; i++) {
            hm.put(groups[i], col[i]);
        }
        try {
            BufferedWriter[] bw = new BufferedWriter[2];
            bw[0] = IOUtils.getTextWriter(outfileS);
            bw[1] = IOUtils.getTextWriter(outfileS1);
            //先写表头
            bw[0].write("TREE_COLORS\nSEPARATOR TAB\nDATA\n");
            bw[1].write("TREE_COLORS\nSEPARATOR TAB\nDATA\n");
            //再写内部的分组，注意文件名字必须和HashMap里的分组名保持一致
            for (int i = 0; i < fs.length; i++) {
                String infileS = fs[i].getAbsolutePath();
                String group = fs[i].getName().split(".txt")[0];
                if(group.equals("NA")) continue;
                BufferedReader br = IOUtils.getTextReader(infileS);
                String temp = null;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    bw[0].write(temp + "\trange\t" + hm.get(group) + "\t" + group);
                    bw[0].newLine();
                    
                    bw[1].write(temp + "\tbranch\t" + hm.get(group) + "\t" + "normal");
                    bw[1].newLine();
                }
                br.close();

            }
            bw[0].flush();
            bw[1].flush();
            bw[0].close();
            bw[1].close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    public void colRangesbyDico() {
        String[] group = {"dicoccum", "dicoccoides"};
        String[] col = {"#76EEC6", "#cd3333"};  // 枫叶红 #cd3333   黄#FFC125  绿#76eec6
        HashMap<String, String> hm = new HashMap<>();
        for (int i = 0; i < group.length; i++) {
            hm.put(group[i], col[i]);
        }
        String infileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/wheatVMapII_AB_S205_germplasmInfo.txt";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/AB/003_tree/002_labelsChange/";
        Arrays.sort(group);
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter[] bw = new BufferedWriter[2];
            //bw[0] = IOUtils.getTextWriter(new File(outfileDirS, "001_addcol_branch_byContient.txt").getAbsolutePath());
            bw[1] = IOUtils.getTextWriter(new File(outfileDirS, "003_colRangebyDico.ABgenome.txt").getAbsolutePath());
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            //bw[0].write("TREE_COLORS\nSEPARATOR TAB\nDATA\n");
            bw[1].write("TREE_COLORS\nSEPARATOR TAB\nDATA\n");
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String DatabaseID = l.get(0);
                String type = l.get(7);
                //bw[0].write(DatabaseID + "\tbranch\t" + hm.get(partContinent));

                int index = Arrays.binarySearch(group, type);
                if (index >= 0) {
                    bw[1].write(DatabaseID + "\trange\t" + hm.get(type) + "\t" + type);
                    bw[1].newLine();
                }
                //bw[0].newLine();
            }
            br.close();
            //bw[0].flush();
            bw[1].flush();
            //bw[0].close();
            bw[1].close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 根据亚洲部洲和其他大洲信息，分组为6组进行着色显示
     */
    public void colRanges() {
        String[] group = {"Oceania", "Africa", "North America", "South America", "Europe", "Asia"};
        //String[] col = {"#F1E1FF","#F4D03F","#F1948A","#5DADE2","#ABEBC6","#239B56","#CD6155","#FF6347","#7B241C","#EBEDEF"};
        String[] col = {"#5DADE2", "#7B241C", "#F1E1FF", "#F4D03F", "#FF9900", "#82C782"};
        HashMap<String, String> hm = new HashMap<>();
        for (int i = 0; i < group.length; i++) {
            hm.put(group[i], col[i]);
        }
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/ABD/001_tree/source/002_merge/All373wheat_ABD_CountryBreedingStatus.txt";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/ABD/001_tree/002_addcolor/";
        Arrays.sort(group);
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter[] bw = new BufferedWriter[2];
            //bw[0] = IOUtils.getTextWriter(new File(outfileDirS, "001_addcol_branch_byContient.txt").getAbsolutePath());
            bw[1] = IOUtils.getTextWriter(new File(outfileDirS, "003_addcol_range_byContient.txt").getAbsolutePath());
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            //bw[0].write("TREE_COLORS\nSEPARATOR TAB\nDATA\n");
            bw[1].write("TREE_COLORS\nSEPARATOR TAB\nDATA\n");
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String DatabaseID = l.get(2);
                String partContinent = l.get(14);
                String continent = l.get(8);
                //bw[0].write(DatabaseID + "\tbranch\t" + hm.get(partContinent));

                int index = Arrays.binarySearch(group, continent);
                if (index >= 0) {
                    bw[1].write(DatabaseID + "\trange\t" + hm.get(continent) + "\t" + partContinent);
                    bw[1].newLine();
                }
                //bw[0].newLine();
            }
            br.close();
            //bw[0].flush();
            bw[1].flush();
            //bw[0].close();
            bw[1].close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 修改树的颜色、标签风格（'normal','bold', 'italic' or 'bold-italic'）、字体大小
     */
    public void TREE_COLORS() {

//        try {
//            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/AB/003_tree/002_labelsChange/002_reaptItem.txt";
//            String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/AB/003_tree/002_labelsChange/003_reaptItemlist.txt";
//            String temp;
//            BufferedReader br = IOUtils.getTextReader(infileS);
//            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
//            while ((temp = br.readLine()) != null) {
//                bw.write(PStringUtils.fastSplit(temp).get(1).substring(0, 5));
//                bw.newLine();
//                bw.write(PStringUtils.fastSplit(temp).get(1).substring(6));
//                bw.newLine();
//            }
//            br.close();
//            bw.flush();
//            bw.close();
//
//        } catch (Exception e) {
//            e.printStackTrace();
//            System.exit(1);
//        }
        try {
            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/AB/003_tree/002_labelsChange/reaptItemlist.txt";
            String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/AB/003_tree/002_labelsChange/002_treeColors.ABgenome.txt";
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("TREE_COLORS\n");
            bw.write("SEPARATOR TAB\n");
            bw.write("DATA\n");
            String temp = null; //read header
            while ((temp = br.readLine()) != null) {
                String databaseID = PStringUtils.fastSplit(temp).get(0);
                //String taxa = PStringUtils.fastSplit(temp).get(4);
                bw.write(databaseID + "\tlabel\t#008000\tbold\t1.5");
                bw.newLine();
                bw.write(databaseID + "\trange\t#D3D3D3\trepeat");
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 修改树node的名字
     */
    public void labels() {
        String infileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/种质信息库/wheatVMapII_AB_S205_germplasmInfo.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/AB/003_tree/002_labelsChange/001_labelschange.ABgenome.txt";
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("LABELS\n");
            bw.write("SEPARATOR TAB\n");
            bw.write("DATA\n");
            String temp = br.readLine(); //read header
            while ((temp = br.readLine()) != null) {
                String databaseID = PStringUtils.fastSplit(temp).get(0);
                String taxa = PStringUtils.fastSplit(temp).get(4);
                bw.write(databaseID + "\t" + databaseID + " " + taxa);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 根据育种状态信息进行 binary标记。 3种标签：Traditional cultivar/Landrace Landrace
     * Advanced/improved cultivar Cultivar Breeding/Research Material
     * 3种颜色：#FFD700 #87CEFA #FF6A6A 3种标记：HV 2 3 状态：hm1 为 landrace的状态标记；如果是
     * Traditional cultivar/Landrace，为1，如果是 Landrace，为0，如果是其他，为-1； hm2 为
     * cultivar的状态标记；如果是 Cultivar，为1，如果是Advanced/improved cultivar，为0；如果是其他，为-1；
     * hm3 为 breeding material的状态标记；如果是 Breeding/Research Material，为1，如果不是，为-1.
     */
    public void binarybyType() {
        String[] type = {"Traditional cultivar/Landrace", "Landrace", "Advanced/improved cultivar", "Cultivar", "Breeding/Research Material", "NA", "Other"};
        HashMap<String, String> hm1 = new HashMap<>();
        HashMap<String, String> hm2 = new HashMap<>();
        HashMap<String, String> hm3 = new HashMap<>();
        //System.out.println(type);
        for (int i = 0; i < type.length; i++) {
            if (type[i].equals("Traditional cultivar/Landrace")) {
                hm1.put(type[i], "1");
                hm2.put(type[i], "-1");
                hm3.put(type[i], "-1");
            }
            if (type[i].equals("Landrace")) {
                hm1.put(type[i], "0");
                hm2.put(type[i], "-1");
                hm3.put(type[i], "-1");
            }
            if (type[i].equals("Cultivar")) {
                hm1.put(type[i], "-1");
                hm2.put(type[i], "1");
                hm3.put(type[i], "-1");
            }
            if (type[i].equals("Advanced/improved cultivar")) {
                hm1.put(type[i], "-1");
                hm2.put(type[i], "0");
                hm3.put(type[i], "-1");
            }
            if (type[i].equals("Breeding/Research Material")) {
                hm1.put(type[i], "-1");
                hm2.put(type[i], "-1");
                hm3.put(type[i], "1");
            }
            if (type[i].equals("Other")) {
                hm1.put(type[i], "-1");
                hm2.put(type[i], "-1");
                hm3.put(type[i], "-1");
            }
            if (type[i].equals("NA")) {
                hm1.put(type[i], "-1");
                hm2.put(type[i], "-1");
                hm3.put(type[i], "-1");
            }
        }

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/ABD/001_tree/source/002_merge/All373wheat_ABD_CountryBreedingStatus.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/ABD/001_tree/002_addcolor/002_binarybyType.txt";
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            bw.write("DATASET_BINARY\n");
            bw.write("SEPARATOR TAB\n");
            bw.write("DATASET_LABEL\tbinary_data\n");
            bw.write("COLOR\t#ff0000\n\n");
            bw.write("FIELD_LABELS\tLandrace\tCultivar\tBreeding/Research Material\n");
            bw.write("FIELD_COLORS\t#FFD700\t#87CEFA\t#FF6A6A\n");
            bw.write("FIELD_SHAPES\t1\t2\t3\n\n");
            bw.write("DATA\n");
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String DatabaseID = l.get(2);
                String status = l.get(13);
                bw.write(DatabaseID + "\t" + hm1.get(status) + "\t" + hm2.get(status) + "\t" + hm3.get(status) + "\n");
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 根据亚洲部洲和其他大洲信息，分组为10组进行着色显示
     */
    public void colbyContinent() {
        String[] group = {"Oceania", "Africa", "North America", "South America", "Europe", "Central Asia", "South Asia", "Western Asia", "East Asia", "NA"};
        //String[] col = {"#F1E1FF","#F4D03F","#F1948A","#5DADE2","#ABEBC6","#239B56","#CD6155","#FF6347","#7B241C","#EBEDEF"};
        String[] col = {"#5DADE2", "#7B241C", "#F1E1FF", "#F4D03F", "#FF9900", "#006600", "#389038", "#82C782", "#CCFFCC", "#EBEDEF"};
        HashMap<String, String> hm = new HashMap<>();
        for (int i = 0; i < group.length; i++) {
            hm.put(group[i], col[i]);
        }
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/ABD/001_tree/source/002_merge/All373wheat_ABD_CountryBreedingStatus.txt";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/ABD/001_tree/002_addcolor";
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter[] bw = new BufferedWriter[2];
            bw[0] = IOUtils.getTextWriter(new File(outfileDirS, "001_addcol_branch_byContient.txt").getAbsolutePath());
            bw[1] = IOUtils.getTextWriter(new File(outfileDirS, "001_addcol_range_byContient.txt").getAbsolutePath());
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            bw[0].write("TREE_COLORS\nSEPARATOR TAB\nDATA\n");
            bw[1].write("TREE_COLORS\nSEPARATOR TAB\nDATA\n");
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String DatabaseID = l.get(2);
                String partContinent = l.get(14);
                bw[0].write(DatabaseID + "\tbranch\t" + hm.get(partContinent) + "\tnormal");
                bw[1].write(DatabaseID + "\trange\t" + hm.get(partContinent) + "\t" + partContinent);
                bw[0].newLine();
                bw[1].newLine();
            }
            br.close();
            bw[0].flush();
            bw[1].flush();
            bw[0].close();
            bw[1].close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 根据育种状态，将373样品分为7类，并赋予其一定颜色： Traditional cultivar/Landrace PaleGreen3	124
     * 205 124	#7CCD7C Advanced/improved cultivar LightGoldenrod1	255 236 139
     * #FFEC8B Breeding/Research Material DarkGreen	0 100 0	#006400 Landrace
     * Green3	0 205 0	#00CD00 Cultivar Goldenrod1	255 193 37	#FFC125 Other
     * LightBlue	173 216 230	#ADD8E6 NA Black	0 0 0	#000000
     */
    public void colbyType() {
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_tree/source/002_merge/All373wheat_ABD_CountryBreedingStatus.txt";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_tree/002_addcolor/";
        BufferedWriter[] bw = new BufferedWriter[2];

        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            bw[0] = IOUtils.getTextWriter(new File(outfileDirS, "001_addcol_branch.txt").getAbsolutePath());
            bw[1] = IOUtils.getTextWriter(new File(outfileDirS, "001_addcol_range.txt").getAbsolutePath());
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String DatabaseID = l.get(2);
                String type = l.get(13);
                if (type.equals("Traditional cultivar/Landrace")) {
                    bw[0].write(DatabaseID + " branch #7CCD7C normal");
                    bw[1].write(DatabaseID + " range #7CCD7C Traditional cultivar/Landrace");
                    bw[0].newLine();
                    bw[1].newLine();
                }
                if (type.equals("Advanced/improved cultivar")) {
                    bw[0].write(DatabaseID + " branch #FFEC8B normal");
                    bw[1].write(DatabaseID + " range #FFEC8B Advanced/improved cultivar");
                    bw[0].newLine();
                    bw[1].newLine();
                }
                if (type.equals("Breeding/Research Material")) {
                    bw[0].write(DatabaseID + " branch #006400 normal");
                    bw[1].write(DatabaseID + " range #006400 Breeding/Research Material");
                    bw[0].newLine();
                    bw[1].newLine();
                }
                if (type.equals("Landrace")) {
                    bw[0].write(DatabaseID + " branch #00CD00 normal");
                    bw[1].write(DatabaseID + " range #00CD00 Landrace");
                    bw[0].newLine();
                    bw[1].newLine();
                }
                if (type.equals("Cultivar")) {
                    bw[0].write(DatabaseID + " branch #FFC125 normal");
                    bw[1].write(DatabaseID + " range #FFC125 Cultivar");
                    bw[0].newLine();
                    bw[1].newLine();
                }
                if (type.equals("Other")) {
                    bw[0].write(DatabaseID + " branch #ADD8E6 normal");
                    bw[1].write(DatabaseID + " range #ADD8E6 Other");
                    bw[0].newLine();
                    bw[1].newLine();
                }
                if (type.equals("NA")) {
                    bw[0].write(DatabaseID + " branch #000000 normal");
                    bw[1].write(DatabaseID + " range #000000 NA");
                    bw[0].newLine();
                    bw[1].newLine();
                }
            }
            br.close();
            bw[0].flush();
            bw[1].flush();
            bw[0].close();
            bw[1].close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
