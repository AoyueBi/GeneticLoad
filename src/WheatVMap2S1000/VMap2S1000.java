package WheatVMap2S1000;

import AoUtils.AoFile;
import AoUtils.CalVCF;
import AoUtils.SplitScript;
import PopulationAnalysis.XPCLR;
import WheatGeneticLoad.FilterVCF2;
import WheatGeneticLoad.VariantsSum;
import daxing.common.IOTool;
import daxing.common.PGF;
import daxing.common.RowTableTool;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import pgl.infra.anno.gene.GeneFeature;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.IntStream;

public class VMap2S1000 {
    public VMap2S1000(){
        /**
         * Finalize VMap 2.0_S1062
         */
//        this.bgzip();
//        this.rename();
//        this.filterSNPtoBi_parallel();

        /**
         * vcf QC
         */
//        new FilterVCF2().QC();
//        new FilterVCF2().statVcfDepth_SD();
//        new FilterVCF2().mkDepthOfVMapII();
//        new FilterVCF2().mkDepthSummary();
//        this.extractVCF(); //提取不同倍性的VCF
//        new FilterVCF2().QC();

//        this.getVCFbyPolyploid();
//        this.getMAFcountfromPop();

        /**
         * gene site annotation
         */
//        this.geneInfo(); //列出所要建立数据库的基因的详细信息表格
        this.snpAnnotationBuild(); //include many methods XXXXXXX
//        new DeleteriousCount();
//        this.filterN_fromVCF();

        /**
         * XPCLR
         */
//        new XPCLR();
//        new DeleteriousXPCLRS1000();




    }

    public void snpAnnotationBuild(){
//        this.mkGeneVCF(); //自己的方法（最终采用）
//        this.mkGeneVCF2(); //来自达兴的方法
//        this.extractInfoFromGeneVCF();
//        this.extractInfoFromGeneVCF_byAoyue();
//        new VariantsSum().mkExonAnnotation2(); //未用
//        new VariantsSum().addAncestral();
//        this.addDAF();
        new VariantsSum().addGerp();
        /**
         * sift 计算
         */
//        new SIFT();
//        this.calFileLine();

//        new VariantsSum().addSift_20210717();
//        new VariantsSum().addDerived_SIFT();
//        new VariantsSum().mergeExonSNPAnnotation();
    }

    /**
     * 获取642倍体中的MAF 个数
     */
    public void getMAFcountfromPop(){

        /**
         * test
         */
//        String infileS = "/Users/Aoyue/Documents/in/chr005.vcf";
//        String outfileS = "/Users/Aoyue/Documents/out/chr005.txt";
//        String taxaInfoDB = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/001_taxaCheck/000_taxaList/VcfIDList_DD.txt";
//        CalVCF.calMAFcountfromPop(infileS,outfileS,taxaInfoDB);

//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/201_VMap2.1";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/039_MAFcount/aabbdd";
//        String logDirS = "/data4/home/aoyue/vmap2/aaPlantGenetics/log_20210810";
//        String taxaInfoDB = "/data4/home/aoyue/vmap2/analysis/038_subsetVCF/VcfIDList_AABBDD.txt";
//        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};


//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/201_VMap2.1";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/039_MAFcount/aabb";
//        String logDirS = "/data4/home/aoyue/vmap2/aaPlantGenetics/log_20210810";
//        String taxaInfoDB = "/data4/home/aoyue/vmap2/analysis/038_subsetVCF/VcfIDList_AABB.txt";
//        String[] chrArr ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};

        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/201_VMap2.1";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/039_MAFcount/dd";
        String logDirS = "/data4/home/aoyue/vmap2/aaPlantGenetics/log_20210810";
        String taxaInfoDB = "/data4/home/aoyue/vmap2/analysis/038_subsetVCF/VcfIDList_DD.txt";
        String[] chrArr ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};

//        for (int i = 0; i < chrArr.length; i++) {
//            String chr = chrArr[i];
//            String infileS = new File(infileDirS,"chr" + chr + "_vmap2.1.vcf").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chr + "_vmap2.1_MAFcount.txt").getAbsolutePath();
//            String logfileS = new File(logDirS,"log_" + new File(outfileS).getName().split(".gz")[0]).getAbsolutePath(); //不管是不是gz结尾，我们只取gz前的部分，妙！
//            System.out.println("java -jar 054_calMAFcountFromPop.jar " +  infileS + " " + outfileS + " " + taxaInfoDB + " > log_20210810_dd.txt");
////            CalVCF.calMAFcountfromPop(infileS,outfileS,taxaInfoDB);
////            System.out.println(chr + "\t finished");
//        }

        SplitScript.splitScript2("/Users/Aoyue/Documents/sh.sh",14,3);
    }




    public void getVCFbyPolyploid(){

        //*** 提取四倍体的VCF ***//
//        String infileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF_Subgenome/500k/chrAB_vamp2.1_500k.vcf.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF_Ploidy/tetra.vcf.gz";
//        List<String> taxaList = AoFile.getStringListwithoutHeader("/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/001_taxaCheck/000_taxaList/VcfIDList_AABB.txt",0);
//        CalVCF.extractVCF(infileS,outfileS,taxaList);

        //*** AABBDD - AB sub ***//
//        String infileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF_Subgenome/500k/chrAB_vmap2.1.500k.vcf.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF_Ploidy/hexaploid_ABsub.vcf.gz";
//        List<String> taxaList = AoFile.getStringListwithoutHeader("/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/001_taxaCheck/000_taxaList/VcfIDList_AABBDD.txt",0);
//        CalVCF.extractVCF(infileS,outfileS,taxaList);

        //*** DD - D sub ***//
        String infileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF_Subgenome/500k/chrD_vmap2.1_500k.vcf.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF_Ploidy/Diploid.vcf.gz";
        List<String> taxaList = AoFile.getStringListwithoutHeader("/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/001_taxaCheck/000_taxaList/VcfIDList_DD.txt",0);
        CalVCF.extractVCF(infileS,outfileS,taxaList);

        //*** AABBDD - D sub ***//
//                String infileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF_Subgenome/500k/chrD_vmap2.1.500k.vcf.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF_Ploidy/hexaploid_Dsub.vcf.gz";
//        List<String> taxaList = AoFile.getStringListwithoutHeader("/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/001_taxaCheck/000_taxaList/VcfIDList_AABBDD.txt",0);
//        CalVCF.extractVCF(infileS,outfileS,taxaList);


        //*** AABBDD ***//
//        String infileDirS = "/data4/home/aoyue/vmap2/daxing/analysis/007_vmap2_1062/001_subsetVmap2.1/001_vmap2_500k/001_byChrID";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/038_subsetVCF/001_hexaploid";
//        String taxaList = "/data4/home/aoyue/vmap2/analysis/038_subsetVCF/VcfIDList_AABBDD.txt";
//        String logDirS = "/data4/home/aoyue/vmap2/aaPlantGenetics/log_049_20210808";
//
//        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
//        for (int i = 0; i < chrArr.length; i++) {
//            String chr = chrArr[i];
//            String infileS = new File(infileDirS,"chr" + chr + "_vmap2.1.500k.vcf.gz").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chr + "_vmap2.1_hexaploid_subset.vcf.gz").getAbsolutePath();
//            String logfileS = new File(logDirS,"log_" + new File(outfileS).getName().split(".gz")[0]).getAbsolutePath(); //不管是不是gz结尾，我们只取gz前的部分，妙！
//            System.out.println("nohup java -jar 049_extractVCF_GL.jar " + infileS + " " + outfileS + " " + taxaList + " > " + logfileS );
//        }

//        SplitScript.splitScript2("/Users/Aoyue/Documents/sh.sh",20,3);




    }


    /**
     * 过滤GeneVCF中alt 含有N的位点
     */
    public void filterN_fromVCF(){
//        String infileDirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/002_geneSNPVCF";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/003_geneSNPVCF_RemoveN";

        String infileDirS = "/data1/publicData/wheat/genotype/VMap/VMap2.0/VMap2.0";
        String outfileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/202_VMap2.0";

        File[] fs = new File(infileDirS).listFiles();
        File[] fs1 = IOUtils.listFilesEndsWith(fs,"vcf.gz");
        List<File> fsList = new ArrayList<>();
        for (int i = 0; i < fs1.length; i++) {
            fsList.add(fs1[i]);
        }
        Collections.sort(fsList);

        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS,f.getName().replaceFirst(".vcf.gz",".vcf")).getAbsolutePath();
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String temp = null;
                int cnttotal = 0;
                int cntsubset = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("##")) {
                        bw.write(temp);
                        bw.newLine();
                    } else if(temp.startsWith("#C")){
                        List<String> l = PStringUtils.fastSplit(temp);
                        bw.write("#CHROM");
                        StringBuilder sb = new StringBuilder();
                        for (int i = 1; i < l.size(); i++) {
                            sb.append("\t").append(l.get(i));
                        }
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    else {
                        cnttotal++;
                        List<String> l = PStringUtils.fastSplit(temp);
                        if (l.get(4).equals("N")) continue; //alt是index为4的列
                        bw.write(temp);
                        bw.newLine();
                        cntsubset++;
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName() + "\t" + cnttotal + "\t" + cntsubset);
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
        // java -jar GeneticLoad.jar > log_removeNfromVCF_20210717.txt 2>&1 &
        // java -jar GeneticLoad.jar > log_removeNfromVMap2.0_20210806.txt 2>&1 &


    }


    public void calFileLine (){
//        String infileDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/002_sift/003_output/001_alt";
//        String infileDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/002_sift/003_output/002_ref";
        String infileDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/006_geneSNPAnnotation_merge";

        File[] fsArray = AoFile.getFileArrayInDir(infileDirS);
        for (int i = 0; i < fsArray.length; i++) {
            String filename = fsArray[i].getName();
            String infileS = fsArray[i].getAbsolutePath();
            String chr = filename.substring(3,6);
            int lineNum =  AoFile.countFileRowNumber_withHeader(infileS);
            System.out.println(chr + "\t" + lineNum);
        }
    }

    public void extractVCF(){

        String infileS = "/Users/Aoyue/project/wheatVMap2_1000/001_germplasm/008_WheatVMap2_GermplasmInfo_20210708.txt";
        String infileS2 = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subset/chrDsubgenome.100k.vcf.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF/diploid.vcf.gz";
        List<String> taxaList = new ArrayList<>();
        RowTable t = new RowTable(infileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            String vcfID = t.getCellAsString(i,0);
            String genomeType = t.getCellAsString(i,4);
            if (genomeType.equals("DD")){
                taxaList.add(vcfID);
            }
        }
        CalVCF.extractVCF(infileS2,outfileS,taxaList);


//        String infileS = "/Users/Aoyue/project/wheatVMap2_1000/001_germplasm/009_WheatVMap2_GermplasmInfo_20210708.txt";
//        String infileS2 = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subset/chrABsubgenome.100k.vcf.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/003_vcfQC/002_subsetVCF/hexaploid.vcf.gz";
//        List<String> taxaList = new ArrayList<>();
//        RowTable t = new RowTable(infileS);
//        for (int i = 0; i < t.getRowNumber(); i++) {
//            String vcfID = t.getCellAsString(i,0);
//            String genomeType = t.getCellAsString(i,4);
//            if (genomeType.equals("AABBDD")){
//                taxaList.add(vcfID);
//            }
//        }
//
//
//        CalVCF.extractVCF(infileS2,outfileS,taxaList);

    }



    /**
     * 老师的方法
     */
    public void addDAF () {
//        String dirS = "/data4/home/aoyue/vmap2/analysis/027_annoDB/002_genicSNP/003_exonSNPAnnotation";
//        String dirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/001_geneSNPByChr";
        String dirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/001_geneSNP_Annotation";
        File[] fs = new File (dirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".txt");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String header = null;
            List<String> recordList = new ArrayList();
            String tem = null;
            try {
                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
                header = br.readLine();
                List<String> l = PStringUtils.fastSplit(header);
                StringBuilder sb = new StringBuilder(header);
                sb.append("\tDAF");
                header = sb.toString();
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    recordList.add(temp);
                }
                br.close();
                BufferedWriter bw = AoFile.writeFile(f.getAbsolutePath());
                bw.write(header);
                bw.newLine();
                float daf = -1;
                String subMajor = null;
                String subMinor = null;
                String ancestral = null;
                String derivedSIFT = null;
                double subMaf = -1;
                for (int i = 0; i < recordList.size(); i++) {
                    sb.setLength(0);
                    sb.append(recordList.get(i)).append("\t");
                    l = PStringUtils.fastSplit(recordList.get(i));
                    ancestral = l.get(9);
                    if (l.get(5).equals(ancestral)) { //major
                        sb.append((float)Double.parseDouble(l.get(7)));
                    }
                    else if (l.get(6).equals(ancestral)) { //minor
                        sb.append((float)(1- Double.parseDouble(l.get(7))));
                    }
                    else sb.append("NA");
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                System.out.println(f.getName());
            }
            catch (Exception e) {
                System.out.println(tem);
                e.printStackTrace();
            }
        });
        // java -Xms50g -Xmx200g -jar PlantGenetics.jar > log_addDAFbyFeisMethod_20200720.txt 2>&1 &
        // java -Xms50g -Xmx200g -jar GeneticLoad.jar > log_addDAFbyFeisMethod_20210716.txt 2>&1 &
        // java -Xms50g -Xmx200g -jar GeneticLoad.jar > log_addDAFbyFeisMethod_20210806.txt 2>&1 &

    }

    /**
     * 由于FEI写的程序中，直接从INFO提取 MAF 和判断 Major Minor 会出现错误，故自己写程序计算该列信息
     */
    public void extractInfoFromGeneVCF_byAoyue () {
        /**
         * inputFile
         */
//        String outDirS = "/Users/Aoyue/Documents/out";
//        String vmapDirS = "/Users/Aoyue/Documents/in";
//        String geneHCFileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/001_geneTable/wheat_v1.1_nonoverlap_addPos.txt.gz";

        String outDirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/001_geneSNP_Annotation";
        String vmapDirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/002_geneSNPVCF";
        String geneHCFileS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/001_geneHC/wheat_v1.1_nonoverlap_addPos.txt.gz";


        /**
         * parameters need to modify
         */
        File[] fs  = AoFile.getFileArrayInDir(vmapDirS);
        List<File> vmapList = Arrays.asList(fs);
        Collections.sort(vmapList);
        AoFile.readheader(geneHCFileS);
        RowTable<String> t = new RowTable<>(geneHCFileS);

        TIntHashSet chrSet = new TIntHashSet(t.getColumnAsIntArray(5)); //get chr的set集合
        int chrNum = chrSet.size(); //chr的个数
        TIntList[] startLists = new TIntList[chrNum]; //list类型的数组，每个数组存放一堆list值
        TIntList[] endLists = new TIntList[chrNum];
        List<String>[] tranLists = new ArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) {
            startLists[i] = new TIntArrayList();
            endLists[i] = new TIntArrayList();
            tranLists[i] = new ArrayList();
        }

        for (int i = 0; i < t.getRowNumber(); i++) {
            int ifunique = Integer.parseInt(t.getCell(i,3));
            if (ifunique == 0) continue; //过滤不是unique的基因
            startLists[Integer.parseInt(t.getCell(i, 5))-1].add(Integer.parseInt(t.getCell(i, 6))); //获取开始位点
            endLists[Integer.parseInt(t.getCell(i, 5))-1].add(Integer.parseInt(t.getCell(i, 7))); //获取终止位点
            tranLists[Integer.parseInt(t.getCell(i, 5))-1].add(t.getCell(i, 4)); //获取基因名字
        }
        vmapList.parallelStream().forEach(f -> {
            int chrIndex = Integer.parseInt(f.getName().substring(3, 6))-1;
            String outfileS = new File (outDirS, f.getName().replaceFirst("_gene_vmap2.1.vcf", "_geneSNPAnno.txt")).getAbsolutePath();
            StringBuilder sb = new StringBuilder();
            sb.append("ID\tChr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tTranscript");
            try {
                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
                BufferedWriter bw = AoFile.writeFile(outfileS);
                bw.write(sb.toString());
                bw.newLine();
                String temp = null;
                List<String> l = null;
                List<String> ll = null;
                List<String> lll = null;
                String info = null;
                int currentPos = -1;
                int posIndex = -1;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#"))continue;
                    sb.setLength(0);
                    l = PStringUtils.fastSplit(temp);
                    currentPos = Integer.parseInt(l.get(1));
                    // 开始判断位点是否在基因中
                    posIndex = startLists[chrIndex].binarySearch(currentPos); //在开始位点搜索是否在
                    if (posIndex < 0) {
                        posIndex = -posIndex-2;
                    }
                    if (posIndex < 0) continue;
                    if (currentPos >= endLists[chrIndex].get(posIndex)) continue;
                    /**
                     * 计算 maf major minor 信息
                     */
                    List<String> lTaxaGeno = new ArrayList<>();
                    String ref = l.get(3);
                    String alt = l.get(4);
                    for (int i = 9; i < l.size(); i++) { //无论有无基因型，都加进去了
                        lTaxaGeno.add(l.get(i));
                    }
                    String[] taxaGenoArray = lTaxaGeno.toArray(new String[lTaxaGeno.size()]);
                    List<String> mafmajorMinorInfo = CalVCF.getMAFMajorMinor(taxaGenoArray,ref,alt);

                    /**
                     * 开始写文件
                     */
                    cnt++;
                    //######## ID ########################### Chr ######################## Pos ####################### Ref
                    sb.append(l.get(2)).append("\t").append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(3));
                    //##################### Alt ##################################### Major ######################################### minor ############################## maf 四位小数
                    sb.append("\t").append(l.get(4)).append("\t").append(mafmajorMinorInfo.get(1)).append("\t").append(mafmajorMinorInfo.get(2)).append("\t").append(mafmajorMinorInfo.get(0));
                    sb.append("\t").append(tranLists[chrIndex].get(posIndex));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getAbsolutePath() + "\t" + cnt  + " is completed.");
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });

        //在HPC上运行： nohup java -jar GeneticLoad.jar > log_extractInfoFromVMap2_20210806.txt 2>&1 &
    }


    public void extractInfoFromGeneVCF () {
        /**
         * inputFile
         */
//        String outDirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/001_geneSNPByChr";
//        String vmapDirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/002_geneSNPVCF";
//        String geneHCFileS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/001_geneHC/wheat_v1.1_nonoverlap_addPos.txt.gz";

        String outDirS = "/Users/Aoyue/Documents/out";
        String vmapDirS = "/Users/Aoyue/Documents/in";
        String geneHCFileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/001_geneTable/wheat_v1.1_nonoverlap_addPos.txt.gz";

        /**
         * parameters need to modify
         */
        int subLength = 200;
        File[] fs  = AoFile.getFileArrayInDir(vmapDirS);
        List<File> vmapList = Arrays.asList(fs);
        Collections.sort(vmapList);
        AoFile.readheader(geneHCFileS);
        RowTable<String> t = new RowTable<>(geneHCFileS);

        TIntHashSet chrSet = new TIntHashSet(t.getColumnAsIntArray(5)); //get chr的set集合
        int chrNum = chrSet.size(); //chr的个数
        TIntList[] startLists = new TIntList[chrNum]; //list类型的数组，每个数组存放一堆list值
        TIntList[] endLists = new TIntList[chrNum];
        List<String>[] tranLists = new ArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) {
            startLists[i] = new TIntArrayList();
            endLists[i] = new TIntArrayList();
            tranLists[i] = new ArrayList();
        }

        for (int i = 0; i < t.getRowNumber(); i++) {
            int ifunique = Integer.parseInt(t.getCell(i,3));
            if (ifunique == 0) continue; //过滤不是unique的基因
            startLists[Integer.parseInt(t.getCell(i, 5))-1].add(Integer.parseInt(t.getCell(i, 6))); //获取开始位点
            endLists[Integer.parseInt(t.getCell(i, 5))-1].add(Integer.parseInt(t.getCell(i, 7))); //获取终止位点
            tranLists[Integer.parseInt(t.getCell(i, 5))-1].add(t.getCell(i, 4)); //获取基因名字
        }
        vmapList.parallelStream().forEach(f -> {
            int chrIndex = Integer.parseInt(f.getName().substring(3, 6))-1;
            String outfileS = new File (outDirS, f.getName().replaceFirst("_vmap2.1.vcf", "_geneSNP.txt")).getAbsolutePath();
            StringBuilder sb = new StringBuilder();
            sb.append("ID\tChr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tTranscript");
            try {
                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
                BufferedWriter bw = AoFile.writeFile(outfileS);
                bw.write(sb.toString());
                bw.newLine();
                String temp = null;
                List<String> l = null;
                List<String> ll = null;
                List<String> lll = null;
                String info = null;
                int currentPos = -1;
                int posIndex = -1;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#"))continue;
                    sb.setLength(0);
                    int currentSub = subLength;
                    if (temp.length() < subLength) { //如果该行的字符长度小于 200，那么最长设置为该行的实际长度
                        currentSub = temp.length();
                    }
                    l = PStringUtils.fastSplit(temp.substring(0, currentSub));
                    currentPos = Integer.parseInt(l.get(1));
                    // 开始判断位点是否在基因中
                    posIndex = startLists[chrIndex].binarySearch(currentPos); //在开始位点搜索是否在
                    if (posIndex < 0) {
                        posIndex = -posIndex-2;
                    }
                    if (posIndex < 0) continue;
                    if (currentPos >= endLists[chrIndex].get(posIndex)) continue;

                    cnt++;
                    //######## ID ########################### Chr ######################## Pos ####################### Ref
                    sb.append(l.get(2)).append("\t").append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(3));
                    //##################### Alt ##################
                    sb.append("\t").append(l.get(4)).append("\t");
                    //start - 弃用提取 INFO 的方法
                    ll = PStringUtils.fastSplit(l.get(7), ";");
                    lll = PStringUtils.fastSplit(ll.get(2).replaceFirst("AD=", ""),",");
                    if (Integer.parseInt(lll.get(0)) > Integer.parseInt(lll.get(1))) {
                        sb.append(l.get(3)).append("\t").append(l.get(4)).append("\t");
                    }
                    else {
                        sb.append(l.get(4)).append("\t").append(l.get(3)).append("\t");
                    }
                    sb.append(ll.get(6).split("=")[1]);
                    //end - 弃用提取 INFO 的方法

                    sb.append("\t").append(tranLists[chrIndex].get(posIndex));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getAbsolutePath() + "\t" + cnt  + " is completed.");
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });

        //在HPC上运行： nohup java -jar GeneticLoad.jar > log_extractInfoFromVMap2_20210716.txt 2>&1 &
    }


    /**
     * 提取高置信度的基因的外显子VCF文件
     */
    public void mkGeneVCF2 () {
        /**
         * inputFile
         */
        String vmapDirS = "/data1/publicData/wheat/genotype/VMap/VMap2.0/VMap2.0"; //modify
        String geneFeatureFileS = "/data1/publicData/wheat/annotation/gene/v1.1/wheat_v1.1_Lulab.pgf"; //modify
        String hcGeneFileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_nonoverlap.txt"; //modify
        String outputDirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/002_geneSNPVCF"; //modify

        RowTableTool<String> rowTable = new RowTableTool<>(hcGeneFileS);
        Predicate<List<String>> p = l->Integer.parseInt(l.get(3))==0;
        rowTable.removeIf(p);

        List<String> geneList = rowTable.getColumn(0);
        PGF pgf = new PGF(geneFeatureFileS);
        Predicate<PGF.Gene> predicate = gene -> !geneList.contains(gene);
        pgf.removeIf(predicate);

        pgf.sortGeneByGeneRange();
        List<File> files = IOTool.getFileListInDirEndsWith(vmapDirS, ".vcf.gz");
        String[] outNames=  files.stream().map(File::getName).map(s -> s.replaceAll("_vmap2.0.vcf.gz","_gene_vmap2.0.vcf")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->{
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bw =IOTool.getWriter(new File(outputDirS, outNames[e]))) {
                String line, subLine;
                int chrID; int pos;
                List<String> temp;
                while ((line= br.readLine()).startsWith("##")){
                    bw.write(line);
                    bw.newLine();
                }
                bw.write(line);
                bw.newLine();
                while ((line=br.readLine())!=null){
                    subLine=line.substring(0, 100);
                    temp=PStringUtils.fastSplit(subLine);
                    chrID =Integer.parseInt(temp.get(0));
                    pos=Integer.parseInt(temp.get(1));
                    int geneIndex=pgf.getGeneIndex(chrID, pos);
                    if (geneIndex < 0) continue;
                    bw.write(line);
                    bw.newLine();
                }
                bw.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
        //java -Xms50g -Xmx200g -jar GeneticLoad.jar > log_mkGeneVCF_20210715.txt 2>&1 &
    }

    /**
     * 提取高置信度的基因的外显子VCF文件
     */
    public void mkGeneVCF () {
        /**
         * inputFile
         */
        String vmapDirS = "/data1/publicData/wheat/genotype/VMap/VMap2.0/VMap2.0"; //modify
        String geneFeatureFileS = "/data1/publicData/wheat/annotation/gene/v1.1/wheat_v1.1_Lulab.pgf"; //modify
        String hcGeneFileS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/001_geneHC/wheat_v1.1_nonoverlap_addPos.txt.gz"; //modify
        String outputDirS = "/data4/home/aoyue/vmap2/analysis/036_annoDB/002_genicSNP/002_geneSNPVCF"; //modify
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        gf.sortGeneByName();

        RowTable<String> t = new RowTable<>(hcGeneFileS);
        TIntHashSet chrSet = new TIntHashSet(t.getColumnAsIntArray(5)); //get chr的set集合
        List<Integer> chrList = new ArrayList<>();
        for (int i = 0; i < chrSet.size(); i++) {
            chrList.add(i+1);
        }
        chrList.parallelStream().forEach(chrID -> {
            String inputVCF = new File(vmapDirS, "chr"+PStringUtils.getNDigitNumber(3, chrID)+"_vmap2.0.vcf.gz").getAbsolutePath();
            String outputVCF = new File (outputDirS, "chr"+PStringUtils.getNDigitNumber(3, chrID)+"_gene_vmap2.1.vcf").getAbsolutePath();
            List<String> geneList = new ArrayList<>();
            List<String> tranList = new ArrayList<>();
            TIntArrayList tranStartList = new TIntArrayList();
            TIntArrayList tranEndList = new TIntArrayList();
            //获取当前 chr 包含的所有 gene 和 trans
            for (int i = 0; i < t.getRowNumber(); i++) {
                int currentChr = Integer.parseInt(t.getCell(i, 5));
                int ifunique = Integer.parseInt(t.getCell(i,3));
                if (ifunique == 0) continue; //过滤不是unique的基因
                if (currentChr < chrID) continue;
                else if (currentChr > chrID) break;
                geneList.add(t.getCell(i, 0));
                tranList.add(t.getCell(i, 4));
                tranStartList.add(t.getCellAsInteger(i,6));
                tranEndList.add(t.getCellAsInteger(i,7));
            }

            //获取每个基因对应的最长转录本的起始终止位点，然后加入大库的 all gene list
            int[] starts = tranStartList.toArray();
            int[] ends = tranEndList.toArray();

            try {
                BufferedReader br = AoFile.readFile(inputVCF);
                BufferedWriter bw = AoFile.writeFile(outputVCF);
                String temp = null;
                while ((temp = br.readLine()).startsWith("##")) {
                    bw.write(temp); bw.newLine();
                }
                bw.write(temp); bw.newLine(); //#CHROM 这一行
                List<String> l = new ArrayList<>();
                int index = -1;
                int pos = -1;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp.substring(0, 100));
                    pos = Integer.parseInt(l.get(1));
                    String alt = l.get(4);
                    if (alt.equals("D") || alt.equals("I") || alt.equals("N")) continue;
                    index = Arrays.binarySearch(starts, pos);
                    if (index < 0) index = -index - 2;
                    if (index < 0) continue;
                    if (pos < ends[index]) {
                        cnt++;
                        bw.write(temp);bw.newLine();
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(chrID + "\t" + cnt  + "\t" + "mkGeneVCF");

            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        //java -Xms50g -Xmx200g -jar GeneticLoad.jar > log_mkGeneVCF_20210715.txt 2>&1 &
        //java -Xms50g -Xmx200g -jar GeneticLoad_2.jar > log_mkGeneVCF_20210806.txt 2>&1 &
    }

    /**
     * 在wheat_v1.1_nonoverlap.txt文件中添加列信息：Chr,TransStart,TransEnd,TranStrand,CDSExonNumber,CDSLength六列信息
     */
    public void geneInfo(){
        String infileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_nonoverlap.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/001_geneTable/wheat_v1.1_nonoverlap_addPos.txt.gz";
        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/wheat/gene/v1.1/wheat_v1.1_Lulab.pgf";
        GeneFeature gf = new GeneFeature (geneFeatureFileS);
        PGF pgf = new PGF(geneFeatureFileS);
        pgf.sortGeneByName();
        gf.sortGeneByName();
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String header = br.readLine();
            bw.write(header + "\tChr\tTransStart\tTransEnd\tTranStrand\tExonNumber\tCDSLength");bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                sb.setLength(0);
                l = PStringUtils.fastSplit(temp);
                cnt++;
                String gene = l.get(0);
                String trans = l.get(4);
                int geneIndex = gf.getGeneIndex(gene);
                int longIndex = gf.getLongestTranscriptIndex(geneIndex);
                String longTrans = gf.getTranscriptName(geneIndex,longIndex);
                if (!trans.equals(longTrans)) System.out.println(temp);
                int chr = gf.getChromosomeOfGene(geneIndex);
                if (chr==0)continue;
                int start = gf.getTranscriptStart(geneIndex,longIndex);
                int end = gf.getTranscriptEnd(geneIndex,longIndex);
                int strand = gf.getTranscriptStrand(geneIndex,longIndex);
                int exonNum = gf.getExonList(geneIndex,longIndex).size();
                int CDSlength = pgf.getCDSLen(geneIndex, longIndex);
                sb.append("\t").append(chr).append("\t").append(start).append("\t").append(end).append("\t").append(strand).append("\t").append(exonNum).append("\t").append(CDSlength);
                bw.write(temp + sb.toString());
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



    public void rename(){
        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
        for (int i = 0; i < chrArr.length; i++) {
            String chr = chrArr[i];
//            System.out.println("mv chr" + chr + ".vmap2.vcf.gz chr" + chr + "_vmap2.0.vcf.gz");
            System.out.println("mv chr" + chr + ".vmap2.vcf.gz.tbi chr" + chr + "_vmap2.0.vcf.gz.tbi");
        }
    }


    /**
     * 构建 VMap2.1 使 alt 保持为1，只含有 A T G C 这几种类型
     */
    public void filterSNPtoBi_parallel() {
        String infileDirS = "/data1/publicData/wheat/genotype/VMap/VMap2.0/VMap2.0";
        String outfileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/201_VMap2.1";
        File[] fs = new File(infileDirS).listFiles();
        File[] fs1 = IOUtils.listFilesEndsWith(fs,"vcf.gz");
        File[] fs2 = IOUtils.listFilesEndsWith(fs,"vcf");
        List<File> fsList = new ArrayList<>();
        for (int i = 0; i < fs1.length; i++) {
            fsList.add(fs1[i]);
        }
        for (int i = 0; i < fs2.length; i++) {
            fsList.add(fs2[i]);
        }

        Collections.sort(fsList);
        System.out.println("Chr\tBiallelicNum\tInsertionNum\tDelectionNum\tTri-alleleNum\tN_Num");
        fsList.parallelStream().forEach(f -> {
            try {
                String chr = f.getName().substring(3, 6); //提取染色体号 001
                int chrint = Integer.parseInt(chr); //将染色体号转化为数字
                String outfileS = new File(outfileDirS, "chr" + chr + "_vmap2.1.vcf").getAbsolutePath();
                BufferedReader br = AoFile.readFile(f.getAbsolutePath());
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String temp = null;
                int biallelicNum = 0;
                int deletion = 0;
                int insertion = 0;
                int triNum = 0;
                int NANum =0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp);
                        bw.newLine();
                    } else {
                        String alt = PStringUtils.fastSplit(temp).get(4);
                        if (alt.contains(",")) {
                            triNum++;
                            continue;
                        }
                        if (alt.equals("D")) {
                            deletion++;
                            continue;
                        }
                        if (alt.equals("I")) {
                            insertion++;
                            continue;
                        }
                        if(alt.equals("N")) {
                            NANum++;
                            continue;
                        }

                        biallelicNum++;
                        bw.write(temp);
                        bw.newLine();
                    }
                }

                br.close();
                bw.flush();
                bw.close();
                System.out.println(chrint + "\t" + biallelicNum + "\t" +  insertion + "\t" + deletion + "\t" + triNum + "\t" + NANum);
//                System.out.println(chrint + "\t" + biallelicNum + "\t" + f.getName() + " is completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }


    /**
     * 对Fei过滤生成的filter2 进行解压缩，用bgzip压缩，再建立索引
     */
    public void bgzip(){
        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
        for (int i = 0; i < chrArr.length; i++) {
            String chr = chrArr[i];
//            System.out.println("nohup gunzip chr" + chr + ".vmap2.vcf.gz 2>&1 &");
//            System.out.println("nohup bgzip -@ 20 chr" + chr + ".vmap2.vcf && tabix -p vcf chr" + chr + ".vmap2.vcf.gz &");
//            System.out.println("nohup bgzip chr" + chr + "_vmap2.0.vcf && tabix -p vcf chr" + chr + "_vmap2.0.vcf.gz &");
            System.out.println("nohup tabix -p vcf chr" + chr + "_vmap2.0.vcf.gz &");

            //测试失败，只能在原来文件夹中压缩
//            System.out.println("nohup bgzip -@ 4 chr" + chr + ".vmap2.vcf -c ../vmap2.0/chr" + chr + ".vmap2.0.vcf.gz && cd ../vmap2.0/ " + "&& tabix -p vcf chr" + chr + ".vmap2.0.vcf.gz &");

        }
    }

}