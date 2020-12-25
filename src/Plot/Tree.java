/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Plot;

import AoUtils.AoFile;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

/**
 *
 * @author Aoyue
 */
public class Tree {

    public Tree() {
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
//        this.colRangebyHexaDiGroup_ABsubgenome();
//        this.removeDot();
//        this.modifyMegaName();
//        this.labels_Asub();
//        this.binarybyContinent();

        /**
         * 处理学博方法的tree
         */

//        this.colRangebyHexaTetraGroup_Asubgenome();
//        this.labels_Asub_xuebo();
//        this.binarybyContinent();

//        this.colRange_gerp();
//        this.textLabel_gerp();
//        this.test();
//        new AoMath().countCaseInGroup("/Users/Aoyue/Documents/test.txt",0);

        this.extractTreeID();




    }

    /**
     * 从氨基酸序列比对文件中提取文件名字
     */
    public void extractTreeID(){
        String infileS = "/Users/Aoyue/Documents/IGDB/04_PHD/未分类笔记/tree/hsp/Fig2a_hsp.txt";
        String outfileS ="/Users/Aoyue/Documents/IGDB/04_PHD/未分类笔记/tree/hsp/Fig2a_hsp_IDlist.txt";
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith(">")) {
                    bw.write(temp.substring(1));
                    bw.newLine();
                    cnt++;
                }
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    public void test(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/006_fromxuebo/002_labels/trash.txt";
        String infile2S = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/006_fromxuebo/002_labels/key.txt";
        HashMap<String, String> hm = new AoFile().getHashMapStringKey(infileS,1,0);
        List<String> l = new AoFile().getStringList(infile2S,0);
        for (int i = 0; i < l.size(); i++) {
            System.out.println(l.get(i) + "\t" + hm.get(l.get(i)));
        }
        String a = hm.get("hhhhh");
        int j = 0;
        System.out.println(a + " test when there is no key");
    }

    /**
     * 将分组信息添加到树的外面
     */
    public void textLabel_gerp(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/006_GERP/source/group";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/006_GERP/002_label/002_textLabelbySubfamily_Asubgenome.txt";

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
        //凤梨科          竹亚科          虎尾草牙科       稻亚科        黍亚科         早熟禾亚科
        String[] groups = {"Bambusoideae", "Bromeliaceae","Chloridoideae","Oryzoideae","Panicoideae","Pooideae"};
        String[] col = { "#00cc00","#0500ec","#33ccff","#ffcc00","#FFA07A","#EE82EE"};
        HashMap<String, String> hm = new HashMap<>();
        for (int i = 0; i < groups.length; i++) {
            hm.put(groups[i], col[i]);
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            //先写表头
            bw.write("DATASET_TEXT\nSEPARATOR TAB\nDATASET_LABEL\texample text dataset\nCOLOR\t#ff0000\nDATA\n");
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
                    bw.write(temp + "\t" + group + "\t-1\t"+ hm.get(group) + "\tnormal\t1\t90"  );
                    bw.newLine();
                }
                br.close();

            }
            bw.flush();
            bw.flush();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 对GERP中的树进行颜色调整
     */
    public void colRange_gerp(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/006_GERP/source/group";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/006_GERP/002_label/001_colRangebySubfamily_Asubgenome.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/006_GERP/002_label/001_colBranchbySubfamily_Asubgenome.txt";

        String reheaderS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";
        RowTable<String> t = new RowTable<>(reheaderS);
        HashMap<String,String> hmtaxa = new HashMap<>();
        for (int i = 0; i < t.getRowNumber() ; i++) {
            String taxa = t.getCell(i,0);
            String taxaID = t.getCell(i,1);
            hmtaxa.put(taxa,taxaID);
        }


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
                            //凤梨科          竹亚科          虎尾草牙科       稻亚科        黍亚科         早熟禾亚科
        String[] groups = {"Bambusoideae", "Bromeliaceae","Chloridoideae","Oryzoideae","Panicoideae","Pooideae"};
        String[] col = { "#00cc00","#0500ec","#33ccff","#ffcc00","#FFA07A","#EE82EE"};
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

                    bw[1].write(temp+ "\tbranch\t" + hm.get(group) + "\t" + "normal");
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
     * 修改树node的名字,先修改4倍体的名字，再修改六倍体的名字
     */
    public void labels_Asub_xuebo() {

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/source/labelChange_Tetraploid.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/006_fromxuebo/002_labels/003_labelcChange.Asub_tetraploid.txt";
        String reheaderS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";
        

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/source/labelChange_Asub_Hexaploid.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/006_fromxuebo/002_labels/003_labelcChange.Asub_hexaploid.txt";
//        String reheaderS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/006_fromxuebo/002_labels/005_labelcChange.Asub.txt";


        RowTable<String> t = new RowTable<>(reheaderS);
        HashMap<String,String> hm = new HashMap<>();
        for (int i = 0; i < t.getRowNumber() ; i++) {
            String taxa = t.getCell(i,0);
            String taxaID = t.getCell(i,1);
            hm.put(taxa,taxaID);
        }

        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("LABELS\n");
            bw.write("SEPARATOR TAB\n");
            bw.write("DATA\n");
            String temp = br.readLine(); //read header
            while ((temp = br.readLine()) != null) {
                String taxa = PStringUtils.fastSplit(temp).get(0);
                String taxaID = hm.get(taxa);
                String country = PStringUtils.fastSplit(temp).get(2);
                String partContinent = PStringUtils.fastSplit(temp).get(3);
                String index = PStringUtils.fastSplit(temp).get(6);
                bw.write(taxa + "\t" + index + " " + taxa + " " +  country + " " + partContinent);
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
     * 根据学博的方法，进行分组，按照倍性添加颜色
     */
    public void colRangebyHexaTetraGroup_Asubgenome() {

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/006_from003";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/006_fromxuebo/002_labels/001_colRangebySubspecies_Asubgenome.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/006_fromxuebo/002_labels/001_colBranchbySubspecies_Asubgenome.txt";

        String reheaderS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";
        RowTable<String> t = new RowTable<>(reheaderS);
        HashMap<String,String> hmtaxa = new HashMap<>();
        for (int i = 0; i < t.getRowNumber() ; i++) {
            String taxa = t.getCell(i,0);
            String taxaID = t.getCell(i,1);
            hmtaxa.put(taxa,taxaID);
        }


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

        String[] groups = {"Ae.tauschii", "Wild_emmer","Domesticated_emmer","Free_threshing_tetraploid","OtherTetraploid","Cultivar","Landrace","OtherHexaploid"};
        String[] col = { "#87cef9","#ffd702","#7f5701","#016699","#00f3ff","#9900ff","#fc6e6e","#fe63c2"};
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



     //根据大洲信息进行 binary标记。
     // 3种颜色：#FFD700 #87CEFA #FF6A6A 3种标记：HV 2 3 状态：
     //hm1 为 landrace的状态标记；如果是Traditional cultivar/Landrace，为1，如果是 Landrace，为0，如果是其他，为-1；
     //hm2 为cultivar的状态标记；如果是 Cultivar，为1，如果是Advanced/improved cultivar，为0；如果是其他，为-1；
     //hm3 为 breeding material的状态标记；如果是 Breeding/Research Material，为1，如果不是，为-1.
    public void binarybyContinent() {

        //建立颜色值的hashMap 亚洲      欧洲       美洲      非洲      大洋洲
        String[] col = {"#82c782","#ff9900","#cf99ff","#7b241c","#5dace2"};
        String[] shape = {"3","1","2","4","6"};
//        String[] type = {"Western Asia", "East Asia","South Asia","Central Asia","Western Europe","Eastern Europe","Southern Europe","Southeast Europe","Central Europe","Northern Europe", "North America","South America","Africa","Oceania"};
        String[] type = {"Asia","Europe","America","Africa","Oceania"};

        for (int i = 0; i < type.length; i++) {
            System.out.println(type[i]);
        }

        HashMap<String,String>  hmcol = new HashMap<>();
        HashMap<String,String>  hmshape = new HashMap<>();
        for (int i = 0; i < type.length; i++) {
            if(type[i].contains("Asia")){
                hmcol.put(type[i],col[0]);
                hmshape.put(type[i],shape[0]);
            }
            if(type[i].contains("Europe")){
                hmcol.put(type[i],col[1]);
                hmshape.put(type[i],shape[1]);
            }
            if(type[i].contains("America")){
                hmcol.put(type[i],col[2]);
                hmshape.put(type[i],shape[2]);
            }
            if(type[i].contains("Africa")){
                hmcol.put(type[i],col[3]);
                hmshape.put(type[i],shape[3]);
            }
            if(type[i].contains("Oceania")){
                hmcol.put(type[i],col[4]);
                hmshape.put(type[i],shape[4]);
            }
        }

        HashMap<String,String>[]  hms = new HashMap[type.length];
        for (int i = 0; i < type.length; i++) { //初始化 hashmap
            hms[i] = new HashMap<>();
        }

        //System.out.println(type);
        for (int i = 0; i < type.length; i++) { //从第一种类型开始循环
            for (int j = 0; j < hms.length; j++) {
                if(type[j].equals(type[i])){ 
                    hms[j].put(type[i], "1");
                }
                else{
                    hms[j].put(type[i], "-1");
                }
            }
            //System.out.println(type[4] + " is " + hms[0].get(type[4]) + " "+ hms[1].get(type[4]) + " " + hms[2].get(type[4]) + " " + hms[3].get(type[4]) + " " + hms[4].get(type[4]) + " " + hms[5].get(type[4]));
        }

//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/003_labels/004_binarybyType.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/003_labels/004_binarybyType_only5continent.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/006_fromxuebo/002_labels/004_binarybyType_only5continent.txt";
        String reheaderS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";

        RowTable<String> t = new RowTable<>(reheaderS);
        HashMap<String,String> hm = new HashMap<>();
        for (int i = 0; i < t.getRowNumber() ; i++) {
            String taxa = t.getCell(i,0);
            String taxaID = t.getCell(i,1);
            hm.put(taxa,taxaID);
        }
        try {
            BufferedReader br = IOUtils.getTextReader(reheaderS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            bw.write("DATASET_BINARY\n");
            bw.write("SEPARATOR TAB\n");
            bw.write("DATASET_LABEL\tbinary_data\n");
            bw.write("COLOR\t#ff0000\n\n");
            bw.write("FIELD_LABELS");
            for (int i = 0; i < type.length; i++) {
                bw.write("\t" + type[i]);
            }
            bw.newLine();
            bw.write("FIELD_COLORS");
            for (int i = 0; i < type.length; i++) {
                bw.write("\t" + hmcol.get(type[i]));

            }
            bw.newLine();
            bw.write("FIELD_SHAPES");
            for (int i = 0; i < type.length; i++) {
                bw.write("\t" + hmshape.get(type[i]));
            }
            bw.newLine();
            bw.write("DATA\n");
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String taxaID = l.get(0);
                String status = l.get(3);
                bw.write(taxaID);
                for (int i = 0; i < hms.length; i++) {
                    String num = hms[i].get(status);
                    bw.write("\t" + num);
                }
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
    public void labels_Asub() {

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/source/labelChange_Tetraploid.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/003_labels/003_labelcChange.Asub_tetraploid.txt";
//        String reheaderS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";
        
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/source/labelChange_Asub_Hexaploid.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/003_labels/003_labelcChange.Asub_hexaploid.txt";
        String reheaderS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";

        RowTable<String> t = new RowTable<>(reheaderS);
        HashMap<String,String> hm = new HashMap<>();
        for (int i = 0; i < t.getRowNumber() ; i++) {
            String taxa = t.getCell(i,0);
            String taxaID = t.getCell(i,1);
            hm.put(taxa,taxaID);
        }

        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("LABELS\n");
            bw.write("SEPARATOR TAB\n");
            bw.write("DATA\n");
            String temp = br.readLine(); //read header
            while ((temp = br.readLine()) != null) {
                String taxa = PStringUtils.fastSplit(temp).get(0);
                String taxaID = hm.get(taxa);
                String country = PStringUtils.fastSplit(temp).get(3);
                String index = PStringUtils.fastSplit(temp).get(5);
                bw.write(taxaID + "\t" + index );
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
     * 目的，将BARLEY中的.....变成空格键，使之能够在mega中使用
     *
     */
    public void removeDot(){
        try {
            String infileS = "/Users/Aoyue/Documents/chrAB.subgenome.txt";
            String outfileS = "/Users/Aoyue/Documents/chrAB.subgenome_RemoveDot.txt";
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            int cnt = 0;
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                for(int i = 0 ; i < temp.length(); i++){
                    if(String.valueOf(temp.charAt(i)).equals(".")){
                        bw.write(" ");
                    }
                    else{
                        bw.write(String.valueOf(temp.charAt(i)));
                    }
                }
                bw.newLine();
                cnt++;
            }
            bw.flush();
            bw.close();
            br.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }
    
    /**
     * 进行A_B subgenome 的分组
     */
    public void colRangebyHexaDiGroup_ABsubgenome() {
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/004_A_Bsubgenome_maf0.01/source/002_subspecies/"; //分组的文件
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/004_A_Bsubgenome_maf0.01/003_labels/001_colRangebyploidy.A_Bsubgenome.txt"; //输出range文件
//        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/004_A_Bsubgenome_maf0.01/003_labels/001_colBranchbyploidy.A_Bsubgenome.txt"; //输出branch文件

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/006_from003";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/003_labels/001_colRangebySubspecies_Asubgenome.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/003_labels/001_colBranchbySubspecies_Asubgenome.txt";

        String reheaderS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";
        RowTable<String> t = new RowTable<>(reheaderS);
        HashMap<String,String> hmtaxa = new HashMap<>();
        for (int i = 0; i < t.getRowNumber() ; i++) {
            String taxa = t.getCell(i,0);
            String taxaID = t.getCell(i,1);
            hmtaxa.put(taxa,taxaID);
        }


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

        String[] groups = {"Ae.tauschii", "Wild_emmer","Domesticated_emmer","Free_threshing_tetraploid","OtherTetraploid","Cultivar","Landrace","OtherHexaploid"};
        String[] col = { "#87cef9","#ffd702","#7f5701","#016699","#00f3ff","#9900ff","#fc6e6e","#fe63c2"};
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
                    bw[0].write(hmtaxa.get(temp) + "\trange\t" + hm.get(group) + "\t" + group);
                    bw[0].newLine();
                    
                    bw[1].write(hmtaxa.get(temp) + "\tbranch\t" + hm.get(group) + "\t" + "normal");
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
//        String infileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/种质信息库/wheatVMapII_AB_S205_germplasmInfo.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/AB/003_tree/002_labelsChange/001_labelschange.ABgenome.txt";
        
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/003_labels/002_labelcChange.Asub.txt";
        
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("LABELS\n");
            bw.write("SEPARATOR TAB\n");
            bw.write("DATA\n");
            String temp = br.readLine(); //read header
            while ((temp = br.readLine()) != null) {
                String databaseID = PStringUtils.fastSplit(temp).get(0);
                String taxaID = PStringUtils.fastSplit(temp).get(1);
//                bw.write(databaseID + "\t" + databaseID + " " + taxaID);
                bw.write(taxaID + "\t" + databaseID);
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
