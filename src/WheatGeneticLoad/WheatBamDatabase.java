/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package WheatGeneticLoad;

import AoUtils.AoFile;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

import AoUtils.AoMath;
import AoUtils.SplitScript;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.IOFileFormat;

/**
 *
 * @author Aoyue
 */
public class WheatBamDatabase {

    public WheatBamDatabase() {
//        this.VMapIAgenome(); //deprecated this method
//        this.VMapIABgenome(); //deprecated this method
//        this.VMapIABDgenome(); //deprecated this method
//        this.VMapIDgenome(); //deprecated this method
        
        //this.VMapIgenomeMethod2();
        //this.mergeVMapIbamDB(); 
        
        //this.wheatJiaoABDgenome();
        //this.renameJiaoABD();
        //this.mergeVMapIJiaobamDB();

        //this.wheatLuABDgenome();
        //this.mergeVMapIJiaoLubamDB();
        //this.renameLuABD();
        
        //this.getCOV();
        //this.getTaxaBamMap();
        
        //this.renameLuS106();
        //this.wheatLu_ABgenome();
        //this.wheatLu_Dgenome();
        //this.mergeVMapI_JiaoABD_LuABD_AB_D_bamDB();
//        this.getTaxaBamMap_Lu_AB_D();

        /**
         * update bam database
         */
//        this.mkmd5();
//        this.checkIfRunOK();
//        SplitScript.splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/009_mkMD5/mkMD5_ABD_bamfile_remaining_20200419.sh",2,72);
//        this.modifyResult();

//        this.mkHashMapfromNEW2OLD();
//        this.getBamTaxaMap();
//        this.getTaxaBamMap_from817bamDB();


        
    }


    /**
     * step1: get taxa set and build taxa list;
     * step2: if bam is for taxa, taxalist add this bam
     */
    public void getTaxaBamMap_from817bamDB(){

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/010_BamHashMap/001_/All817_WGS_WheatBamDatabase_20200420.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/010_BamHashMap/001_/003_taxa2BamsMap.txt";
        String[] taxaSetArray = AoFile.getStringArraybySet(infileS,1);
        Arrays.sort(taxaSetArray);
        List<String>[] taxaBams = new List[taxaSetArray.length];
        for (int i = 0; i < taxaBams.length; i++) {
            taxaBams[i] = new ArrayList<>();
        }

        try{
            BufferedReader br = AoFile.readFile(infileS);
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            while((temp=br.readLine()) != null){
                l = PStringUtils.fastSplit(temp);
                String bam = l.get(3);
                String taxa = l.get(1);
                int index = Arrays.binarySearch(taxaSetArray,taxa);
                if (index <0){
                    System.out.println(temp);
                }
                taxaBams[index].add(bam);
            }
            br.close();

            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxa\t\"Bams(A list of bams of the taxon, seperated by the delimiter of Tab)\"");
            bw.newLine();
            for (int i = 0; i < taxaBams.length; i++) {
                String taxa = taxaSetArray[i];
                bw.write(taxa);
                for (int j = 0; j < taxaBams[i].size(); j++) {
                    String bam = taxaBams[i].get(j);
                    bw.write("\t" + bam);
                }
                bw.newLine();
            }
            bw.flush(); bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * step1: get all bams name and build hashmap(key:bam value:taxa)
     * step2: get value from vmap2 sheet except AA genome, and then get value from WheatBamDatabase on AA genome
     * step3: complete bamdatabase (DataBaseID	Taxa	Accessions	Genome-type	Bam-Path	Insert-size(bp)	Sequencing-platform	Coverage	DataSource)
     * NAFU: insert size: approximately 500 bp. read length was 150 bp, Illumina HiSeq X Ten platform; coverage ~ 8×
     * VMapII: read length was 100 bp
     * Jiao: read length was 150 bp
     * VMapI: read length was 150bp
     */
    public void getBamTaxaMap(){
        String vmap2fileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/001_toFeiLu/wheatVMapII_germplasmInfo_20191225.txt";
        String bamDBfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/008_step8_mergeVMapI_JiaoABD_LuABD_AB_D/All725WheatBamDatabase_20190716.txt";
        String bamS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/010_BamHashMap/001_/001_new2old_NameHashMap.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/010_BamHashMap/001_/002_newBam2TaxaMap.txt";

        System.out.println("****************************************");
        RowTable<String> t = new RowTable<>(bamS);
        List<String> oldBamList = t.getColumn(1);
        Collections.sort(oldBamList);
        System.out.println("Totally " + oldBamList.size() + " bams in /data3/wgs/bam.");

        HashMap<String,String> hmBamTaxa = new HashMap<String,String>();
        //先判断 vmap2中的bam是否都在 oldBamList中, 将在 oldBamList中的都建立HashMap
        t = new RowTable<>(vmap2fileS);
        List<String> bamListinVMap2 = new ArrayList<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            String bam = t.getCell(i,0);
            String taxa = t.getCell(i,4);
            int index = Collections.binarySearch(oldBamList,bam);
            if (index<0){
                System.out.println(bam + " " + taxa);
                continue;
            }
            hmBamTaxa.put(bam,taxa);
            bamListinVMap2.add(bam);
        }
        Collections.sort(bamListinVMap2);
        System.out.println(hmBamTaxa.size() + " bam taxa map has been got from vmap2");

        //获取不在vmap2中的bamTaxa HashMap
        t = new RowTable<>(bamDBfileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            String bam = t.getCell(i,0);
            String taxa = t.getCell(i,1);
            int index = Collections.binarySearch(bamListinVMap2,bam);
            if (index<0){
                System.out.println(bam + " " + taxa);
                hmBamTaxa.put(bam,taxa);
            }
        }
        System.out.println(hmBamTaxa.size() + " bam taxa map has been got from vmap2 and bamDB");

        //开始写出文件
        List<String> nafuTaxaList = AoFile.getStringList("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/010_BamHashMap/001_/NAFU_test.txt",0);
        HashMap<String,String> hmOld2newbam = AoFile.getHashMapStringKey(bamS,1,0);
        try {
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write("BamDataBaseID\tTaxa\tGenome-type\tBam-Path\tInsert-size(bp)\tRead-length\tSequencing-platform\tCoverage\tDataSource");
            bw.newLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            OldBamDB oldbamDB = new OldBamDB();
            int insertSize = 500;
            int readlength = 150;
            String sequencingPlatform = "Illumina HiSeq X Ten";
            String cov = "8X";
        //* NAFU: insert size: approximately 500 bp. read length was 150 bp, Illumina HiSeq X Ten platform; coverage ~ 8×

            String genometype = null;
            String path = null;
            for (int i = 0; i < oldBamList.size(); i++) {
                String oldBam = oldBamList.get(i);
                String newBam = hmOld2newbam.get(oldBam);
                String taxa = hmBamTaxa.get(oldBam);

                if (newBam.startsWith("A_")){
                    genometype = "AA";
                    path = "/data3/wgs/bam/A/" + newBam+ ".bam";
                }
                if (newBam.startsWith("AB_")){
                    genometype = "AABB";
                    path = "/data3/wgs/bam/AB/"+ newBam+ ".bam";
                }
                if (newBam.startsWith("ABD_") || newBam.startsWith("CS_") ){
                    genometype = "AABBDD";
                    path = "/data3/wgs/bam/ABD/"+ newBam+ ".bam";
                }
                if (newBam.startsWith("D_")){
                    genometype = "DD";
                    path = "/data3/wgs/bam/D/"+ newBam + ".bam";
                }
//bw.write("BamDataBaseID\tTaxa\tGenome-type\tBam-Path\tInsert-size(bp)\tRead-length\tSequencing-platform\tCoverage\tDataSource");
                int index = Collections.binarySearch(nafuTaxaList,taxa);
                if (index > -1){ //说明是属于西农的材料 CS不在其中
                    bw.write(newBam + "\t" + taxa + "\t" + genometype + "\t" + path + "\t" + insertSize + "\t" + readlength + "\t" + sequencingPlatform + "\t" + cov + "\tNorthwest Agriculture and Forestry University");
                    bw.newLine();
                }
                else{ //说明不属于 NAFU的材料
                    if (newBam.startsWith("CS_mp")){
                        bw.write(newBam + "\tCS" + "\t" + genometype + "\t" + path + "\t2-4kb, 5-7kb and 8-10kb" + "\t160" + "\tIllumina HiSeq 2500" + "\t8X" + "\tSRR5893663");
                        bw.newLine();
                    }
                    else if (newBam.startsWith("CS_sg")){
                        bw.write(newBam + "\tCS_3X" + "\t" + genometype + "\t" + path + "\t500" + "\t100 or 150" + "\tIllumina Genome Analyzer IIx" + "\t3X" + "\tERP003210");
                        bw.newLine();
                    }else{
                        bw.write(newBam + "\t" + taxa + "\t" + genometype + "\t" + path + "\t" + oldbamDB.getInsertsize(oldBam) + "\t" + oldbamDB.getReadlength(oldBam) + "\t" + oldbamDB.getSequencingPlatform(oldBam) + "\t" + oldbamDB.getCoverage(oldBam) + "\t" + oldbamDB.getDataSource(oldBam));
                        bw.newLine();
                    }
                }
            }

            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    class OldBamDB{
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/008_step8_mergeVMapI_JiaoABD_LuABD_AB_D/All725WheatBamDatabase_20190716.txt";
        private HashMap<String,Integer> hmBamInsertsize = new HashMap<>();
        private HashMap<String,Integer> hmBamReadlength = new HashMap<>();
        private HashMap<String,String> hmBamSequencingPlatform = new HashMap<>();
        private HashMap<String,String> hmBamCoverage = new HashMap<>();
        private HashMap<String,String> hmBamDataSource = new HashMap<>();



        public OldBamDB(){
            this.readBamFile();
        }

        public void readBamFile(){
            RowTable<String> t = new RowTable<>(infileS);
            int readlength = Integer.MIN_VALUE;
            for (int i = 0; i < t.getRowNumber(); i++) {
                String bam = t.getCell(i,0);
                int insertSize = t.getCellAsInteger(i,5);
                String sequencingPlatform = t.getCell(i,6);
                String coverage = t.getCell(i,7);
                String dataSource = t.getCell(i,8);

                hmBamInsertsize.put(bam,insertSize);
                hmBamSequencingPlatform.put(bam,sequencingPlatform);
                hmBamCoverage.put(bam,coverage);
                hmBamDataSource.put(bam,dataSource);
                if (sequencingPlatform.equals("BGISEQ500")){
                    readlength = 100;
                }
                else{
                    readlength = 150;
                }
                hmBamReadlength.put(bam,readlength);
            }
        }

        public int getInsertsize(String bam){
           return hmBamInsertsize.get(bam);
        }

        public int getReadlength(String bam){
            return hmBamReadlength.get(bam);
        }
        public String getSequencingPlatform(String bam){
            return hmBamSequencingPlatform.get(bam);
        }

        public String getCoverage(String bam){
            return hmBamCoverage.get(bam);
        }

        public String getDataSource(String bam){
            return hmBamDataSource.get(bam);
        }
    }

    /**
     * 将计算的md5结果，路径去掉，保持 MD5value  taxa.bam 格式
     */
    public void modifyResult(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/009_mkMD5/out_bamfile_md5";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/009_mkMD5/out";
        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS, f.getName()).getAbsolutePath();
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String temp = null;
                List<String> l = new ArrayList<>();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp,"  ");
//                    String name = l.get(1).substring(19); //should modify should modify should modify should modify
                    String name = l.get(1).substring(17);
                    System.out.println(name);
                    bw.write(l.get(0) + "  " + name);
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
//                System.out.println(f.getName() + "\tis completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
            }
        });

    }

    /**
     * ABD基因组计算md5的程序被停掉，查看哪些没有跑完，把脚本写出来
     */
    public void checkIfRunOK(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/009_mkMD5/mkMD5_ABD_bamfile_20200418.sh";
        String md5DirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/009_mkMD5/out_bamfile_md5";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/009_mkMD5/mkMD5_ABD_bamfile_remaining_20200419.sh";
        File[] fs = AoFile.getFileArrayInDir(md5DirS);
        List<String> nameList = new ArrayList<>();
        for (int i = 0; i < fs.length; i++) {
            nameList.add(fs[i].getName());
        }
        Collections.sort(nameList);

        try{
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp;
            while((temp = br.readLine()) != null){
                String name = PStringUtils.fastSplit(temp," ").get(3).substring(34);
                System.out.println(name);
                int index = Collections.binarySearch(nameList,name);
                if (index>-1)continue; //搜到了就不写出来了
                bw.write(temp);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }


    }


    /**
     * 对bam文件库进行md5的建立
     */
    public void mkmd5(){
//        for (int i = 1; i < 92; i++) {
//            String num = PStringUtils.getNDigitNumber(4,i);
//            String name = "A_" + num + ".bam";
//            System.out.println("md5sum " + name + " > " + name + ".md5" );
//        }

//        for (int i = 1; i < 232; i++) {
//            String num = PStringUtils.getNDigitNumber(4,i);
//            String name = "AB_" + num + ".bam";
//            System.out.println("md5sum " + name + " > " + name + ".md5" );
//        }

//        String infileDirS = "/data3/wgs/bam/ABD";
//        String outfileDirS = "/data1/home/aoyue/out_bamfile_md5";
//        for (int i = 1; i < 435; i++) {
//            String num = PStringUtils.getNDigitNumber(4,i);
//            String name = "ABD_" + num + ".bam";
//            String infileS = new File(infileDirS,name).getAbsolutePath();
//            String outfileS = new File(outfileDirS,name + ".md5").getAbsolutePath();
//            System.out.println("md5sum " + infileS + " > " + outfileS );
//        }
//        String cs1infileS = new File(infileDirS,"CS_mp_2018_8X.bam").getAbsolutePath();
//        String cs1outfileS = new File(outfileDirS,"CS_mp_2018_8X.bam.md5").getAbsolutePath();
//        String cs2infileS = new File(infileDirS,"CS_sg_2014_3X.bam").getAbsolutePath();
//        String cs2outfileS = new File(outfileDirS,"CS_sg_2014_3X.bam.md5").getAbsolutePath();
//
//        System.out.println("md5sum " + cs1infileS + " > " + cs1outfileS );
//        System.out.println("md5sum " + cs2infileS + " > " + cs2outfileS );

//        SplitScript.splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/009_mkMD5/mkMD5_ABD_bamfile_20200418.sh",10,44);

        /******************************* D bam ********************************************/
//        String infileDirS = "/data3/wgs/bam/D";
//        String outfileDirS = "/data1/home/aoyue/out_bamfile_md5";
//        for (int i = 1; i < 61; i++) {
//            String num = PStringUtils.getNDigitNumber(4,i);
//            String name = "D_" + num + ".bam";
//            String infileS = new File(infileDirS,name).getAbsolutePath();
//            String outfileS = new File(outfileDirS,name + ".md5").getAbsolutePath();
//            System.out.println("md5sum " + infileS + " > " + outfileS );
//        }

//        SplitScript.splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/009_mkMD5/D_mkmd5.sh",10,6);
    }

    public void mkHashMapfromNEW2OLD(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/010_BamHashMap/bam_Old2New.map.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/010_BamHashMap/001_/001_new2old_NameHashMap.txt";
        RowTable<String> t = new RowTable<>(infileS);
        HashMap<String,String> hm = new HashMap<>();
        try {
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write("NewName\tOldName");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) {
                String oldFileS = t.getCell(i,0);
                String newFileS = t.getCell(i,1);
                if (oldFileS.endsWith("bai"))continue;
                // A A0001.rmdup.bam
                // AB B0115.rmdup.bam    NAFU S13.markdup.bam
                //ABD TW0059.rmdup.bam   NAFU ALS.sort.dedup.bam
                // D D0030.rmdup.bam     NAFU A2.markdup.bam
                // new name is endwith ".bam"
                String oldname = new File(oldFileS).getName();
                String newname = new File(newFileS).getName();
                if (oldname.endsWith(".rmdup.bam")){
                    oldname = oldname.split(".rmdup.bam")[0];
                }
                if (oldname.endsWith(".markdup.bam")){
                    oldname = oldname.split(".markdup.bam")[0];
                }
                if (oldname.endsWith(".sort.dedup.bam")){
                    oldname = oldname.split(".sort.dedup.bam")[0];
                }
                if(oldname.equals("CSmarkdup.bam")){
                    oldname = "CS_8X";
                }
                if(oldname.equals("CS.rmdup.bam")){
                    oldname = "CS_8X";
                }
                newname = newname.split(".bam")[0];
                hm.put(oldname,newname);
                System.out.println(oldname + "\t" + newname);
                bw.write(newname + "\t" + oldname);
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println(hm.size());

        }catch (Exception e){
            e.printStackTrace();
            System.out.print(1);
        }

    }
    
    
    public void getTaxaBamMap_Lu_AB_D(){
        /**********************VMapI_AB******************************/
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/003_fastCall/000_uploadParameters/001_calCOV/calCOV_Lu_AB.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/003_fastCall/000_uploadParameters/002_getTaxaBamMap/taxaBamMap_Lu_AB.txt";
//        String bampathDirS = "/data2/aoyue/vmap2_AB_project/003_mergefile/AB/";
        
        /**********************VMapI_D******************************/
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/003_fastCall/000_uploadParameters/001_calCOV/calCOV_Lu_D.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/003_fastCall/000_uploadParameters/002_getTaxaBamMap/taxaBamMap_Lu_D.txt";
        String bampathDirS = "/data1/home/aoyue/vmap2_AB_project/003_mergefile/D/";

        
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxa\tCoverage-Of-All-Bams\t\"Bams(A list of bams of the taxon, seperated by the delimiter of Tab)\"");
            bw.newLine();
            String temp = br.readLine();
            while((temp=br.readLine()) != null){
                String databaseID = PStringUtils.fastSplit(temp).get(0);
                String cov = PStringUtils.fastSplit(temp).get(1);
                bw.write(databaseID + "\t" + cov + "\t" + bampathDirS + databaseID + ".rmdup.bam");
                bw.newLine();
            }
            br.close(); bw.flush(); bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * 将所有文件合并成一个数据库
     */
    public void mergeVMapI_JiaoABD_LuABD_AB_D_bamDB(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/001_mergeVMapI/VMapI_allBam.txt";
        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/002_step2_jiao/Jiao_ABD_bam.txt";
        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/004_step4_luABD/Lu_ABD_bam.txt";
        String infileS3 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/006_step6_Lu_AB/Lu_AB_bam.txt";
        String infileS4 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/007_step7_Lu_D/Lu_D_bam.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/008_step8_mergeVMapI_JiaoABD_LuABD_AB_D/All725WheatBamDatabase_temp.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/008_step8_mergeVMapI_JiaoABD_LuABD_AB_D/All725WheatBamDatabase.txt";
        try{
            //只读入表头
            BufferedReader br = IOUtils.getTextReader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1_method2/VMapI_A_bam.txt");
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(br.readLine());bw.newLine();
            br.close();
            br = IOUtils.getTextReader(infileS);
            String temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                bw.write(temp); bw.newLine();
            }
            br.close();

            br=IOUtils.getTextReader(infileS1);
            temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                bw.write(temp); bw.newLine();
            }
            br.close();
            
            br=IOUtils.getTextReader(infileS2);
            temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                bw.write(temp); bw.newLine();
            }
            
            br=IOUtils.getTextReader(infileS3);
            temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                bw.write(temp); bw.newLine();
            }
            
            br=IOUtils.getTextReader(infileS4);
            temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                bw.write(temp); bw.newLine();
            }
            
            br.close();bw.flush();bw.close();
            
            RowTable t = new RowTable(outfileS);
            t.sortAsText(0);
            t.writeTextTable(outfileS1, IOFileFormat.Text);
        
        new File(outfileS).delete();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    
    
    public void wheatLu_Dgenome(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/Lu_VMapII_D/Lu_D_SeqIDAccession.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/007_step7_Lu_D/Lu_D_bam_temp.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/007_step7_Lu_D/Lu_D_bam.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("DataBaseID\tTaxa\tAccessions\tGenome-type\tBam-Path\tInsert-size(bp)\tSequencing-platform\tCoverage\tDataSource"); bw.newLine();
            String temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                String seqID = PStringUtils.fastSplit(temp).get(0);
                String taxa = PStringUtils.fastSplit(temp).get(3);
                String accessions = PStringUtils.fastSplit(temp).get(1);
                String databaseID = PStringUtils.fastSplit(temp).get(2);
                StringBuilder sb = new StringBuilder();
                sb.append(databaseID).append("\t").append(taxa).append("\t").append(accessions).append("\t").append("D\t").append("/data3/wgs/bam/D/" + databaseID + ".rmdup.bam").append("\t").
                        append("300\t").append("BGISEQ500\t").append("10X\tLuLab");
                bw.write(sb.toString()); bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        RowTable t = new RowTable(outfileS);
        t.sortAsText(0);
        t.writeTextTable(outfileS1, IOFileFormat.Text);
        
        new File(outfileS).delete();
    }
    
    public void wheatLu_ABgenome(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/Lu_VMapII_AB/Lu_AB_SeqIDAccession.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/006_step6_Lu_AB/Lu_AB_bam_temp.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/006_step6_Lu_AB/Lu_AB_bam.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("DataBaseID\tTaxa\tAccessions\tGenome-type\tBam-Path\tInsert-size(bp)\tSequencing-platform\tCoverage\tDataSource"); bw.newLine();
            String temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                String seqID = PStringUtils.fastSplit(temp).get(0);
                String taxa = PStringUtils.fastSplit(temp).get(3);
                String accessions = PStringUtils.fastSplit(temp).get(1);
                String databaseID = PStringUtils.fastSplit(temp).get(2);
                StringBuilder sb = new StringBuilder();
                sb.append(databaseID).append("\t").append(taxa).append("\t").append(accessions).append("\t").append("AB\t").append("/data3/wgs/bam/AB/" + databaseID + ".rmdup.bam").append("\t").
                        append("300\t").append("BGISEQ500\t").append("10X\tLuLab");
                bw.write(sb.toString()); bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        RowTable t = new RowTable(outfileS);
        t.sortAsText(0);
        t.writeTextTable(outfileS1, IOFileFormat.Text);
        
        new File(outfileS).delete();
    }
    
    public void renameLuS106(){
        //String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/Lu_VMapII_AB/Lu_AB_SeqIDAccession.txt";
        //String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/006_step6_Lu_AB/renameLu_AB.sh";
        
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/Lu_VMapII_D/Lu_D_SeqIDAccession.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/007_step7_Lu_D/renameLu_Dmd5.sh";
        /**
         * SeqID	Accessions	DatabaseID	Taxa
            AT08079	CItr 3686	B0126	CItr3686
            AT08080	CItr 14822	B0127	CItr14822
            AT08081	CItr 14824	B0128	CItr14824
         */
        try{
            ///data3/wgs/bam/ABD
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            
            String temp = br.readLine(); // header
            while((temp = br.readLine()) != null){
                String seqID = PStringUtils.fastSplit(temp).get(0);
                String databaseID = PStringUtils.fastSplit(temp).get(2);
                StringBuilder sb = new StringBuilder();
//                sb.append("mv ").append(seqID).append(".rmdup.bam ").append(databaseID).append(".rmdup.bam").
//                        append(" \n").append("mv ").append(seqID).append(".rmdup.bam.bai ").append(databaseID).append(".rmdup.bam.bai");

                sb.append("mv ").append(seqID).append(".rmdup.bam.md5.txt ").append(databaseID).append(".rmdup.bam.md5.txt");
                        
                bw.write(sb.toString()); bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * Get taxaBamMap from clean.fq.gz data
     * 
     * === goal file ===
     * Taxa	Coverage-Of-All-Bams	"Bams(A list of bams of the taxon, seperated by the delimiter of Tab)"
     * TW0001	3	/data3/wgs/bam/ABD/TW0001.rmdup.bam
     * 
     * === central file ===
     * 以Lu_ABD_SeqIDAccession.txt为核心如下，根据DatabaseID找SeqID,和Bam-Path，建立2个HashMap
     * SeqID	Accessions	DatabaseID	Taxa
     * BT01373	CItr 1517	TW0174	CItr1517
     * 
     * ==== cov file ===
     * SeqID	Cov
     * BT01339	11.49206897
     * 
     * === db file ===
     * DataBaseID	Taxa	Accessions	Genome-type	Bam-Path	Insert-size(bp)	Sequencing-platform	Coverage	DataSource
     * A0001	TRI18435	TRI 18435	A	/data3/wgs/bam/A/A0001.rmdup.bam	350	NovaSeq 6000	3X	LuLab
     */
    public void getTaxaBamMap(){
        
        /**********************LuLab_ABD******************************/
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/Lu_ABD/Lu_ABD_SeqIDAccession.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/003_fastCall/000_uploadParameters/002_getTaxaBamMap/taxaBamMap_Lu_ABD.txt";
//        String covfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/003_fastCall/000_uploadParameters/001_calCOV/calCOV_Lu_ABD.txt";
//        String dbfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/005_step5_mergeVMapIandJiaoandLu/All619WheatBamDatabase_20190531.txt";
        
        /**********************JiaoLab_ABD******************************/
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/Jiao_ABD/Jiao_ABD_SeqIDAccession.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/003_fastCall/000_uploadParameters/002_getTaxaBamMap/taxaBamMap_Jiao_ABD.txt";
//        String covfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/003_fastCall/000_uploadParameters/001_calCOV/calCOV_Jiao_ABD.txt";
//        String dbfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/005_step5_mergeVMapIandJiaoandLu/All619WheatBamDatabase_20190531.txt";
        
        /**********************VMapI_ABD******************************/
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/VMapI_ABD.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/003_fastCall/000_uploadParameters/002_getTaxaBamMap/taxaBamMap_VMapI_ABD.txt";
        String covfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/003_fastCall/000_uploadParameters/001_calCOV/calCOV_VMapI_ABD.txt";
        String dbfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/005_step5_mergeVMapIandJiaoandLu/All619WheatBamDatabase_20190531.txt";
        
        /**********************VMapI_AB******************************/
        
        
        
        /**********************VMapI_D******************************/
        
        
        
        HashMap<String,Double> hmSeqCov = new HashMap<>();
        HashMap<String,String> hmDataBaseIDBampath = new HashMap<>();
        
        RowTable<String> t = new RowTable<>(covfileS);
        for(int i=0;i<t.getRowNumber();i++){
            hmSeqCov.put(t.getCell(i, 0), t.getCellAsDouble(i, 1));
        }
        
        t= new RowTable(dbfileS);
        for(int i =0; i<t.getRowNumber();i++){
            hmDataBaseIDBampath.put(t.getCell(i,0), t.getCell(i, 4));
        }
        
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxa\tCoverage-Of-All-Bams\t\"Bams(A list of bams of the taxon, seperated by the delimiter of Tab)\"");
            bw.newLine();
            String temp = br.readLine();
            while((temp=br.readLine()) != null){
                String seqID = PStringUtils.fastSplit(temp).get(0);
                String databaseID = PStringUtils.fastSplit(temp).get(2);
                bw.write(databaseID + "\t" + hmSeqCov.get(seqID) + "\t" + hmDataBaseIDBampath.get(databaseID));
                bw.newLine();
            }
            br.close(); bw.flush(); bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * Get cov from rmdup.bam file
     */
    
    public void getCOV (){
        String infileDirS = "/Users/Aoyue/Documents/cov";
        String outfileS = "/Users/Aoyue/Documents/stat_cov2.txt";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".COV.txt");
        Arrays.sort(fs);
        try{
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("DataBaseID\tCoverage-Of-All-Bams\n");
            for(int i = 0; i < fs.length; i++){
                String infileS = fs[i].getAbsolutePath();
                BufferedReader br = IOUtils.getTextReader(infileS);
                String temp ;
                List<Integer> l = new ArrayList<>();
                HashMap<Integer,String> hm = new HashMap<>();
                while((temp = br.readLine()) != null){
                    String cov = PStringUtils.fastSplit(temp).get(1);
                    int fre = Integer.parseInt(PStringUtils.fastSplit(temp).get(2));
                    hm.put(fre,cov);
                    l.add(fre);
                }
               Collections.sort(l,Collections.reverseOrder());
               String dataBaseID = fs[i].getName().split("_1_10M20M")[0];
               int max = Collections.max(l); String cov = hm.get(max);
               int covN = Integer.parseInt(cov);
//               if( covN == 1){
//                   cov = hm.get(l.get(1));
//               }
               bw.write(dataBaseID + "\t" + cov + "\n");
               br.close();
            }
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void renameLuABD(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/Lu_ABD/Lu_ABD_SeqIDAccession.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/004_step4_luABD/renameLu.sh";
        /**
         * SeqID	Accessions	DatabaseID	Taxa
            BT01373	CItr 1517	TW0174	CItr1517
            BT01374	PI 46041	TW0175	PI46041
            BT01375	PI 61693	TW0176	PI61693
         */
        try{
            ///data3/wgs/bam/ABD
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            
            String temp = br.readLine(); // header
            while((temp = br.readLine()) != null){
                String seqID = PStringUtils.fastSplit(temp).get(0);
                String databaseID = PStringUtils.fastSplit(temp).get(2);
                StringBuilder sb = new StringBuilder();
                sb.append("mv /data3/wgs/wheat200bam/").append(seqID).append("/").append(seqID).append(".sort.dedup.bam /data3/wgs/bam/ABD/").append(databaseID).append(".rmdup.bam")
                        .append(" && ").
                        append("mv /data3/wgs/wheat200bam/").append(seqID).append("/").append(seqID).append(".sort.dedup.bai /data3/wgs/bam/ABD/").append(databaseID).append(".rmdup.bam.bai");
                bw.write(sb.toString()); bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void mergeVMapIJiaoLubamDB(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/001_mergeVMapI/VMapI_allBam.txt";
        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/002_step2_jiao/Jiao_ABD_bam.txt";
        String infileS2 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/004_step4_luABD/Lu_ABD_bam.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/005_step5_mergeVMapIandJiaoandLu/All619WheatBamDatabase.txt";
        try{
            //只读入表头
            BufferedReader br = IOUtils.getTextReader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1_method2/VMapI_A_bam.txt");
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(br.readLine());bw.newLine();
            br.close();
            br = IOUtils.getTextReader(infileS);
            String temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                bw.write(temp); bw.newLine();
            }
            br.close();

            br=IOUtils.getTextReader(infileS1);
            temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                bw.write(temp); bw.newLine();
            }
            br.close();
            
            br=IOUtils.getTextReader(infileS2);
            temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                bw.write(temp); bw.newLine();
            }
            
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void wheatLuABDgenome(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/Lu_ABD/Lu_ABD_SeqIDAccession.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/004_step4_luABD/Lu_ABD_bam_temp.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/004_step4_luABD/Lu_ABD_bam.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("DataBaseID\tTaxa\tAccessions\tGenome-type\tBam-Path\tInsert-size(bp)\tSequencing-platform\tCoverage\tDataSource"); bw.newLine();
            String temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                String seqID = PStringUtils.fastSplit(temp).get(0);
                String taxa = PStringUtils.fastSplit(temp).get(3);
                String accessions = PStringUtils.fastSplit(temp).get(1);
                String databaseID = PStringUtils.fastSplit(temp).get(2);
                StringBuilder sb = new StringBuilder();
                sb.append(databaseID).append("\t").append(taxa).append("\t").append(accessions).append("\t").append("ABD\t").append("/data3/wgs/bam/ABD/" + databaseID + ".rmdup.bam").append("\t").
                        append("300\t").append("BGISEQ500\t").append("10X\tLuLab");
                bw.write(sb.toString()); bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        RowTable t = new RowTable(outfileS);
        t.sortAsText(0);
        t.writeTextTable(outfileS1, IOFileFormat.Text);
        
        new File(outfileS).delete();
    }
    
    public void mergeVMapIJiaobamDB(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/001_mergeVMapI/VMapI_allBam.txt";
        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/002_step2_jiao/Jiao_ABD_bam.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/003_step3_VMapIandJiao/AllWheatBamDatabase.txt";
        try{
            //只读入表头
            BufferedReader br = IOUtils.getTextReader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1_method2/VMapI_A_bam.txt");
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(br.readLine());bw.newLine();
            br.close();
            br = IOUtils.getTextReader(infileS);
            String temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                bw.write(temp); bw.newLine();
            }
            br.close();

            br=IOUtils.getTextReader(infileS1);
            temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                bw.write(temp); bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void renameJiaoABD(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/Jiao_ABD/Jiao_ABD_SeqIDAccession.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/002_step2_jiao/renameJiao.sh";
        /**
         * SeqID	Accessions	DatabaseID
            AFG-L1	PI 80741	TW0055
            AFG-L2	PI 127098	TW0056
            AFG-L3	PI 366569	TW0057
         */
        try{
            ///data3/wgs/bam/ABD
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            
            String temp = br.readLine(); // header
            while((temp = br.readLine()) != null){
                String seqID = PStringUtils.fastSplit(temp).get(0);
                String databaseID = PStringUtils.fastSplit(temp).get(2);
                StringBuilder sb = new StringBuilder();
                sb.append("mv /data3/wgs/bam/ABD/").append(seqID).append(".rmdup.bam /data3/wgs/bam/ABD/").append(databaseID).append(".rmdup.bam").
                        append(" && ").append("mv /data3/wgs/bam/ABD/").append(seqID).append(".rmdup.bam.bai /data3/wgs/bam/ABD/").append(databaseID).append(".rmdup.bam.bai");
                bw.write(sb.toString()); bw.newLine();
            }
            br.close();bw.flush();bw.close();
            
            
            
            
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    

    /**
     * 修改焦老师的文件名，并且进行数据库的建立，再与VMapI进行合并；
     */
    public void wheatJiaoABDgenome(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/Jiao_ABD/Jiao_ABD_SeqIDAccession.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/002_step2_jiao/Jiao_ABD_bam_temp.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/002_step2_jiao/Jiao_ABD_bam.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("DataBaseID\tTaxa\tAccessions\tGenome-type\tBam-Path\tInsert-size(bp)\tSequencing-platform\tCoverage\tDataSource"); bw.newLine();
            String temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                String seqID = PStringUtils.fastSplit(temp).get(0);
                String taxa = PStringUtils.fastSplit(temp).get(3);
                String accessions = PStringUtils.fastSplit(temp).get(1);
                String databaseID = PStringUtils.fastSplit(temp).get(2);
                StringBuilder sb = new StringBuilder();
                sb.append(databaseID).append("\t").append(taxa).append("\t").append(accessions).append("\t").append("ABD\t").append("/data3/wgs/bam/ABD/" + databaseID + ".rmdup.bam").append("\t").
                        append("350\t").append("Illumina HiSeq X\t").append("10X\tJiaoLab");
                bw.write(sb.toString()); bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        RowTable t = new RowTable(outfileS);
        t.sortAsText(0);
        t.writeTextTable(outfileS1, IOFileFormat.Text);
        
        new File(outfileS).delete();
    }
    
    /**
     * 通过同时读取VMapI_A.txt VMapI_AB.txt VMapI_ABD.txt VMapI_D.txt 4个文件，进行数据库的生成，一步到位。
     */
    public void VMapIgenomeMethod2(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1_method2/";
        
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesStartsWith(fs, "VMapI_");
        for (int i=0; i<fs.length;i++){
            System.out.println(fs[i].getName());
        }
        System.out.println(fs.length);
        
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(file -> {
            String genome = file.getName().split("I_")[1].split(".t")[0];
            System.out.println(genome + "   genome is being executed");
            String infileS = file.getAbsolutePath();
            String outfileS = new File(outfileDirS,file.getName().split(".txt")[0] + "_bam_temp.txt" ).getAbsolutePath();
            String outfileS1 = new File(outfileDirS,file.getName().split(".txt")[0] + "_bam.txt" ).getAbsolutePath();
            try{
                BufferedReader br = IOUtils.getTextReader(infileS);
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write("DataBaseID\tTaxa\tAccessions\tGenome-type\tBam-Path\tInsert-size(bp)\tSequencing-platform\tCoverage\tDataSource"); bw.newLine();
                String temp = br.readLine(); //表头
                while((temp = br.readLine()) != null){
                    String taxa = PStringUtils.fastSplit(temp).get(2);
                    String seqID = PStringUtils.fastSplit(temp).get(0);
                    String accessions = PStringUtils.fastSplit(temp).get(1);
                    String databaseID = null;
                    if(seqID.startsWith("B0")){
                        databaseID = seqID.replaceFirst("B0", "B00");
                    }
                    else if (seqID.startsWith("B1")) {
                        databaseID = seqID.replaceFirst("B1", "B01");
                    }
                    else if(seqID.startsWith("A0")){
                        databaseID = seqID.replaceFirst("A0", "A00");
                    }
                    else if(seqID.startsWith("TW0")){
                        databaseID = seqID.replaceFirst("TW0", "TW00");
                    }
                    else if(seqID.startsWith("D0")){
                        databaseID = seqID.replaceFirst("D0", "D00");
                    }
                    StringBuilder sb = new StringBuilder();
                    sb.append(databaseID).append("\t").append(taxa).append("\t").append(accessions).append("\t").append(genome).append("\t").append("/data3/wgs/bam/" + genome + "/" + 
                            databaseID + ".rmdup.bam").append("\t").
                            append("350\t").append("NovaSeq 6000\t").append("3X\tLuLab");
                    bw.write(sb.toString()); bw.newLine();
                }
                br.close();bw.flush();bw.close();
            }
            catch(Exception e){
                e.printStackTrace();
                System.exit(1);
            }
            
            RowTable t = new RowTable(outfileS);
            t.sortAsText(0);
            t.writeTextTable(outfileS1, IOFileFormat.Text);

            new File(outfileS).delete();
        });
        
    }
    
    /**
     * 将生成的A AB ABD 和 D bam文件数据库合并起来，成为一个文件。
     * 1.建立表头，写入头文件；
     * 2.建立数组，将文件名字按顺序读入；
     * 3.进行for循环，合并文件。
     */
    public void mergeVMapIbamDB(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1_method2";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/001_mergeVMapI/VMapI_allBam.txt";
        String[] genomes = {"A", "AB", "ABD", "D"};
        try{
            //只读入表头
            BufferedReader br = IOUtils.getTextReader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1/VMapI_A_bam.txt");
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(br.readLine());bw.newLine();
            br.close();
            //按文件顺序合并
            for(int i =0; i < genomes.length; i++){
                String infileS = new File(infileDirS,"VMapI_" + genomes[i] + "_bam.txt").getAbsolutePath();
                br = IOUtils.getTextReader(infileS);
                String temp = br.readLine(); //表头
                while((temp = br.readLine()) != null){
                    bw.write(temp); bw.newLine();
                }
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
    }
    
    public void VMapIDgenome(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/VMapI_D.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1/VMapI_D_bam_temp.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1/VMapI_D_bam.txt";

        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("DataBaseID\tSeqID\tAccessions\tGenome-type\tBam-Path\tInsert-size(bp)\tSequencing-platform\tCoverage\tDataSource"); bw.newLine();
            String temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                String seqID = PStringUtils.fastSplit(temp).get(0);
                String accessions = PStringUtils.fastSplit(temp).get(1);
                String databaseID = seqID.replaceFirst("D0", "D00");
                StringBuilder sb = new StringBuilder();
                sb.append(databaseID).append("\t").append(seqID).append("\t").append(accessions).append("\t").append("D\t").append("/data3/wgs/bam/D/").append("\t").
                        append("350\t").append("NovaSeq 6000\t").append("3X\tVMapI");
                bw.write(sb.toString()); bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        RowTable t = new RowTable(outfileS);
        t.sortAsText(0);
        t.writeTextTable(outfileS1, IOFileFormat.Text);
        
        new File(outfileS).delete();
    }
    
    public void VMapIABDgenome(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/VMapI_ABD.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1/VMapI_ABD_bam_temp.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1/VMapI_ABD_bam.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("DataBaseID\tSeqID\tAccessions\tGenome-type\tBam-Path\tInsert-size(bp)\tSequencing-platform\tCoverage\tDataSource"); bw.newLine();
            String temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                String seqID = PStringUtils.fastSplit(temp).get(0);
                String accessions = PStringUtils.fastSplit(temp).get(1);
                String databaseID = seqID.replaceFirst("TW0", "TW00");
                StringBuilder sb = new StringBuilder();
                sb.append(databaseID).append("\t").append(seqID).append("\t").append(accessions).append("\t").append("ABD\t").append("/data3/wgs/bam/ABD/").append("\t").
                        append("350\t").append("NovaSeq 6000\t").append("3X\tVMapI");
                bw.write(sb.toString()); bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        RowTable t = new RowTable(outfileS);
        t.sortAsText(0);
        t.writeTextTable(outfileS1, IOFileFormat.Text);
        
        new File(outfileS).delete();
    }
    
    public void VMapIABgenome(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/VMapI_AB.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1/VMapI_AB_bam_temp.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1/VMapI_AB_bam.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("DataBaseID\tSeqID\tAccessions\tGenome-type\tBam-Path\tInsert-size(bp)\tSequencing-platform\tCoverage\tDataSource"); bw.newLine();
            String temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                String seqID = PStringUtils.fastSplit(temp).get(0);
                String accessions = PStringUtils.fastSplit(temp).get(1);
                String databaseID;
                if(seqID.startsWith("B0")){
                    databaseID = seqID.replaceFirst("B0", "B00");
                }
                else{
                    databaseID = seqID.replaceFirst("B1", "B01");
                }
                StringBuilder sb = new StringBuilder();
                sb.append(databaseID).append("\t").append(seqID).append("\t").append(accessions).append("\t").append("AB\t").append("/data3/wgs/bam/AB/").append("\t").
                        append("350\t").append("NovaSeq 6000\t").append("3X\tVMapI");
                bw.write(sb.toString()); bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        RowTable t = new RowTable(outfileS);
        t.sortAsText(0);
        t.writeTextTable(outfileS1, IOFileFormat.Text);
        
        new File(outfileS).delete();
    }
    
    /**
     * 1.get the name of each SeqID, and then rename it!
     * 2.Make table like DataBaseID	SeqID	Accessions	Genome-type	Bam-Path	Insert-size	Sequencing-platform	Coverage
     * 
     * 将文本读进去，建立HashMap
     */
    public void VMapIAgenome(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/VMapI_A.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1/VMapI_A_bam_temp.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1/VMapI_A_bam.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("DataBaseID\tSeqID\tAccessions\tGenome-type\tBam-Path\tInsert-size(bp)\tSequencing-platform\tCoverage\tDataSource"); bw.newLine();
            String temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                String seqID = PStringUtils.fastSplit(temp).get(0);
                String accessions = PStringUtils.fastSplit(temp).get(1);
                String databaseID = seqID.replaceFirst("A0", "A00");
                StringBuilder sb = new StringBuilder();
                sb.append(databaseID).append("\t").append(seqID).append("\t").append(accessions).append("\t").append("A\t").append("/data3/wgs/bam/A/").append("\t").
                        append("350\t").append("NovaSeq 6000\t").append("3X\tVMapI");
                bw.write(sb.toString()); bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        RowTable t = new RowTable(outfileS);
        t.sortAsText(0);
        t.writeTextTable(outfileS1, IOFileFormat.Text);
        
        new File(outfileS).delete();
    }
}
