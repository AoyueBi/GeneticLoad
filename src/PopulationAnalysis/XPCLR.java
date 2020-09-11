package PopulationAnalysis;

import AoUtils.*;

import analysis.wheat.VMap2.VMapDBUtils;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;

import pgl.infra.table.ColumnTable;
import pgl.infra.table.RowTable;
import pgl.infra.utils.Dyad;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import smile.stat.Stat;
import tech.tablesaw.api.IntColumn;
import tech.tablesaw.api.Table;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;
import pgl.graph.tSaw.TablesawUtils;

public class XPCLR {
    public XPCLR(){
//        this.convertCoordinate();
//        this.test();
//        this.addgeneticPos();
//        this.calDensity("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/002_snp/chr001.subgenome.maf0.01.SNP_bi.subset.pos.Base.txt.gz",1,3,2000000,2000000,"/Users/Aoyue/Documents/test.txt");

//        this.getGenotypeXPCLR();
//        this.script_getGenotypeXPCLR();
        /**
         * 多线程打包测试 进行input数据的输入
         */
//        this.getGenotypeXPCLR_parallele();
//        this.getGenotypeXPCLR_parallele_tetra();
//        this.getGenotypeXPCLR_parallele_diploid();
        /**
         * 最终执行脚本文件
         */
//        this.script_XPCLR();
//        this.script_XPCLR_tetraploid();
//        new SplitScript().splitScript2("/Users/Aoyue/Documents/xpclr_20200301.sh",12,4);
//        new SplitScript().splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/006_script/003_WEDE/xpclr_tetraploid_WE_DE_20200310.sh",14,2); //14,2

//        SplitScript.splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/006_script/004_FTTDE/xpclr_tetraploid_FTT_DE_20200312.sh",14,2);
        /**
         * snp density 统计
         */
//        this.script_calSNPdensity();
//        this.mergeTxt();
//        new SplitScript().splitScript2("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/006_script/sh_xpclr_hexaploid20200224.sh",14,3);

//        this.statisticSNPdensity();

        /**
         * 结果处理：合并，转换坐标，
         */
//        this.mergeTxt2();
//        this.convertXPCLRCoordinate(); //for manhatton plot && change chr pos


//        this.convertXPCLRCoordinate2(); //将结果不进行坐标转换，只添加表头，把信息不完全的行删除,进行topK做准备
//        this.sortbyXPCLR(); //已和上一步连起来运行
//        this.getTopK(); //已和上一步连起来运行
//        this.getSelectedPos(); //采用这种方法来获取受选择区域的基因列表

//        this.addGeneID(); //已放弃
//        this.addGeneID_onlyGridPos(); //已放弃


//        this.checkTopGeneDistribution();


        /**
         * **************************************** VMap2.0 new version calculation ********************************************************
         */

//        this.checkAnnotationDBisinExonVCF(); //确定annotation的位点都在exonVCF中
//        this.getExonVCFbyPloidy();

//        this.mkSNPfile_hexaploid(); //分很多步骤
//        this.mkSNPfile_tetraploid();
        this.getXPCLRscript("abd");
//        this.getXPCLRscript("ab");
        //对exon位点数进行计数
//        CountSites.countSites_singleStream("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/002_snp_file");

//        this.X(); //对XPCLR结果进行初处理
//        this.window();

//        this.test();

//        this.getAlleleCount();



    }

    public void getAlleleCount(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/000_exonVCF/chr013_exon_vmap2.1.vcf.gz";
        String objectfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/000_pop_bySubspecies/Cultivar.txt";
        String reffileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/000_pop_bySubspecies/Landrace_EU.txt";
        String outfileS = "/Users/Aoyue/Documents/chr013_exon_vmap2.1.alleleCount.txt";
        this.getAlleleCountXPCLR(infileS,objectfileS,reffileS,outfileS);
    }


    public void getAlleleCountXPCLR(String infileS, String objectfileS, String reffileS, String outfileS){
        List<String> queryTaxal = new AoFile().getStringListwithoutHeader(objectfileS,0); //获取pop列表
        List<String> query2Taxal = new AoFile().getStringListwithoutHeader(reffileS,0); //获取pop列表

        try {
            BufferedReader br = new AoFile().readFile(infileS);
            BufferedWriter bw = new AoFile().writeFile(outfileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            List<Integer> indexObject = new ArrayList<>();
            List<Integer> indexRef = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##"))continue;
                if (temp.startsWith("#C")){
                    l = PStringUtils.fastSplit(temp);
                    for (int i = 9; i < l.size(); i++) {
                        String taxon = l.get(i);

                        int index = Collections.binarySearch(queryTaxal, taxon);
                        if (index > -1) { //当找到列表中的taxa时，写列表中的taxa信息
                            indexObject.add(i);
                        }

                        int index2 = Collections.binarySearch(query2Taxal, taxon);
                        if (index2 > -1) { //当找到列表中的taxa时，写列表中的taxa信息
                            indexRef.add(i);
                        }
                    }
                    Collections.sort(indexObject);
                    Collections.sort(indexRef);

//                    System.out.println("Finish find the pop index from vcffile.");
                }
                if (!temp.startsWith("#")) {
                    l = PStringUtils.fastSplit(temp);
                    List<String> lGeno = new ArrayList<>();
                    List<String> lGeno2 = new ArrayList<>();

                    for (int i = 0; i < indexObject.size(); i++) { //无论有无基因型，都加进去了
                        lGeno.add(l.get(indexObject.get(i)));
                    }


                    for (int i = 0; i < indexRef.size(); i++) { //无论有无基因型，都加进去了
                        lGeno2.add(l.get(indexRef.get(i)));
                    }

                    String[] GenoArray = lGeno.toArray(new String[lGeno.size()]);
                    String[] Geno2Array = lGeno2.toArray(new String[lGeno2.size()]);

                    String geno = this.getAlleleCountInfo(GenoArray);
                    String geno2 = this.getAlleleCountInfo(Geno2Array);

                    bw.write(geno + "\t" + geno2);
                    bw.newLine();
                    cnt++;
                }
            }
            System.out.println(cnt + " SNP " + new File(infileS).getName() + " is completed at " + outfileS);
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * object Pop: 0 0 1 1 0 0 1 1
     * ref Pop: 0 0 0 0 1 0
     * Then, object include 4 individuals, total allele count in object is 8, allele 1 count is 4
     *          ref include 3 individuals, total allele count in object is 6, allele 1 count is 1
     *
     * object Pop: 0 0 9 9 0 0 1 1
     * ref Pop: 0 0 0 0 0 0
     * Then, object include 4 individuals, total allele count in object is 6, allele 1 count is 2
     *        ref include 3 individuals, total allele count in object is 6, allele 1 count is 0
     *
     * @param genoArray
     * @return
     */
    public String getAlleleCountInfo(String[] genoArray) {
        String out = null;
        int totalCounts = 0;
        int allele1Count = 0;


        List<String> tempList = null;
        List<String> temList = null;
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < genoArray.length; i++) {
            tempList = PStringUtils.fastSplit(genoArray[i], ":"); //tempList是包含基因型GT AD PL的集合
            temList = PStringUtils.fastSplit(tempList.get(0), "/"); //temList是包含基因型拆分后的集合
            String ref1 = temList.get(0);
            String alt2 = temList.get(1);
            if (ref1.equals(".")){
                ref1 = "9";
                alt2 = "9";
            }else {
                totalCounts++;totalCounts++;
                if (ref1.equals("1")){
                    allele1Count++;
                }
                if (alt2.equals("1")){
                    allele1Count++;
                }
            }
        }

        sb.append(totalCounts).append("\t").append(allele1Count);
        out = sb.toString();
        return out;
    }



    /**
     * 从XPCLR的结果中获取每条染色体的完成度，及最后一行到第几M
     */
    public void test(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/006_output/102_0.0001_100_500";
        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.stream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                BufferedReader br = AoFile.readFile(infileS);
                String header = br.readLine();
                String temp = null;
                String pos = null;
                List<String> l = new ArrayList<>();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    pos = l.get(3);
                }
                String chr = f.getName().substring(3,6);
                System.out.println(chr + "\t" + pos);
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    /**
     * 将标准化的结果进行滑窗处理
     */
    public void window(){

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/006_output/Cultivar_VS_Landrace_EU_exonRegion_0.0001_100_500.xpclr.txt.gz";
        AoFile.readheader(infileS);
        int chrColumn = 8;
        int posIndex = 9;
        int valueIndex = 10;
        double window = 100000;
        double step = 50000;
        String name = new File(infileS).getName().split(".txt")[0] + "_" + window + "window_" + step + "step.txt.gz";
        String parent = new File(infileS).getParent();
        String outfileS = new File(parent,name).getAbsolutePath();

        System.out.println("nohup java -jar 051_AoWindowScan.jar " + infileS + " " + chrColumn + " " + posIndex + " " + valueIndex + " " + window + " " + step + " " + outfileS + " &" );
    }

    /**
     * 将原始结果添加参考基因组的坐标，并对xpclr标准化
     */
    public void X(){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/006_output/002_0.0001_100_500";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/006_output/001_0.0001_100_100000";
//        String outfileDirS = AoString.autoOutfileDirS(infileDirS);

        String infileDirS = "/Users/Aoyue/Documents/001";
        String outfileDirS = AoString.autoOutfileDirS(infileDirS);

        new File(outfileDirS).mkdirs();
        System.out.println(outfileDirS);
        int a=3;
        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
//            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/006_output/chr001_Clutivar_VS_Landrace_EU_exonRegion_0.0001_100_100000.xpclr.txt";
//            String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/008_xpclrHighlevel/ztest/chr001_Clutivar_VS_Landrace_EU_exonRegion_0.0001_100_100000.xpclr.txt";

            String infileS = f.getAbsolutePath();
            String outfileS = new File(outfileDirS,f.getName()).getAbsolutePath();

            List<String> recordList = new ArrayList<>();
            TDoubleArrayList valueList = new TDoubleArrayList();
            try {
                BufferedReader br = new AoFile().readFile(infileS);
                String temp = null;
                List<String> l = new ArrayList<>();
                int cnt = 0;

                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp," ");
                    if (l.size() < 7)continue; //如果该行的list小于7，说明没有读写完全就终止了，是最后一行，应该去除
                    cnt++;
                    //2 217 1152 21716707.000000 0.651300 584.597428 0.000000
                    double v = Double.parseDouble(l.get(5));
                    valueList.add(v);
                    recordList.add(temp);
                }
                br.close();

                //开始处理标准化
                double[] valueArray = valueList.toArray(new double[valueList.size()]);
                double[] normalizedScore = AoMath.NormalizeScore(valueArray);

                BufferedWriter bw = new AoFile().writeFile(outfileS);
                bw.write("CHROM\tGrid\tN_SNPs\tPOS\tGenetic_pos\tXPCLR_score\tMax_s\tChrRef\tPosRef\tXPCLR_100score");
                bw.newLine();
                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < recordList.size(); i++) {
                    temp = recordList.get(i);
                    l = PStringUtils.fastSplit(temp," ");
                    String chrS = l.get(0);
                    String posS = l.get(3);
                    int chrID = Integer.parseInt(chrS);
                    int posID = (int) Double.parseDouble(posS);
                    String chrref = RefV1Utils.getChromosome(chrID,posID);
                    int posref = RefV1Utils.getPosOnChromosome(chrID,posID);
                    sb.setLength(0);
                    sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(2)).append("\t").
                            append(l.get(3)).append("\t").append(l.get(4)).append("\t").append(l.get(5)).
                            append("\t").append(l.get(6)).append("\t").append(chrref).append("\t").append(posref).append("\t").append(normalizedScore[i]);

                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                System.out.println(cnt + "\t" + new File(infileS).getName() + " is completed.");
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

        });

        String outfileS = new File(new File(infileDirS).getParent(),fsList.get(0).getName().substring(7)).getAbsolutePath();
        AoFile.mergeTxt(outfileDirS,outfileS);
    }


    /**
     *  准备六倍体的输入文件
     */
    public void mkSNPfile_hexaploid(){
//        String pop1fileS = "";
//        String pop2fileS = "";
//        String infileDirS = "";
//        String outParentS = "";

//        String pop1fileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/000_pop_bySubspecies/Cultivar.txt";//goal
//        String pop2fileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/000_pop_bySubspecies/Landrace_EU.txt"; //ref

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/000_exonVCF";
//        String outParentS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid";


        /////// 挑取少量个体，测试chr013是否能跑通
        String pop1fileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/000_pop_bySubspecies_smallPop/Cultivar.txt";
        String pop2fileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/000_pop_bySubspecies_smallPop/Landrace_EU.txt";
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/000_exonVCF";
        String outParentS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/test";

        //根据输出文件1，子目录输出
        String outfileDirS = new File(outParentS,"001_chrPosRefAlt").getAbsolutePath(); new File(outfileDirS).mkdirs();
        String outfileDirS2 = new File(outParentS,"002_snp_file").getAbsolutePath(); new File(outfileDirS2).mkdirs();
        String outfileDirS3 = new File(outParentS,"004_calDensity").getAbsolutePath(); new File(outfileDirS3).mkdirs();
        String outfileDirS4 = new File(outParentS,"005_merged003").getAbsolutePath(); new File(outfileDirS4).mkdirs();
        String outfileDirS5 = new File(outParentS,"003_genotype").getAbsolutePath();new File(outfileDirS5).mkdirs();


//        this.step1_mkSNPfile(infileDirS,outfileDirS);
//        this.step2_mkSNPfile(outfileDirS,outfileDirS2);

        this.getGenotype_parallele(infileDirS,pop1fileS,outfileDirS5);
        this.getGenotype_parallele(infileDirS,pop2fileS,outfileDirS5);

        //        this.getSNPdensity_hexaploid(outfileDirS,outfileDirS3);
//        this.mergeTXT(outfileDirS3,outfileDirS4);

//        this.calDensity("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/001_chrPosRefAlt/chr042_exon_vmap2.1.pos.Base.txt.gz",1,2,100000,100000,"/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/001_chrPosRefAlt/chr042.txt");
    }



    /**
     *
     */
    public void mkSNPfile_tetraploid(){
        String pop1fileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/000_pop_bySubspecies/Domesticated_emmer.txt";//goal
        String pop2fileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/000_pop_bySubspecies/Wild_emmer.txt"; //ref

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/005_tetraploid/000_exonVCF";
        String outParentS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/005_tetraploid";

        //根据输出文件1，子目录输出
        String outfileDirS = new File(outParentS,"001_chrPosRefAlt").getAbsolutePath(); new File(outfileDirS).mkdirs();
        String outfileDirS2 = new File(outParentS,"002_snp_file").getAbsolutePath(); new File(outfileDirS2).mkdirs();
        String outfileDirS3 = new File(outParentS,"003_calDensity").getAbsolutePath(); new File(outfileDirS3).mkdirs();
        String outfileDirS4 = new File(outParentS,"004_merged003").getAbsolutePath(); new File(outfileDirS4).mkdirs();
        String outfileDirS5 = new File(outParentS,"005_genotype").getAbsolutePath();new File(outfileDirS5).mkdirs();


//        this.step1_mkSNPfile(infileDirS,outfileDirS);
//        this.step2_mkSNPfile(outfileDirS,outfileDirS2);
//        this.getSNPdensity_tetraploid(outfileDirS,outfileDirS3);
//        this.mergeTXT(outfileDirS3,outfileDirS4);
//        this.getGenotype_parallele(infileDirS,pop1fileS,outfileDirS5);
//        this.getGenotype_parallele(infileDirS,pop2fileS,outfileDirS5);
        this.getXPCLRscript("ab");

    }




    public void getXPCLRscript(String ploidy){

        String infileDirS = null;
        String outfileDirS = null;
        String snpInfoDirS = null;
        String logDirS = null;
        String[] chrArr;

        String pop1 = null;
        String pop2 = null;
        String gwin = null;
        String snpWin = null;
        String gridSize = null;

        if (ploidy.equals("abd")){
            String parentDirS = "/data4/home/aoyue/vmap2/analysis/030_XPCLR/004_hexaploid/test";
            infileDirS = new File(parentDirS,"003_genotype").getAbsolutePath();
            outfileDirS = new File(parentDirS,"006_output").getAbsolutePath();
            snpInfoDirS = new File(parentDirS,"002_snp_file").getAbsolutePath();
            logDirS = new File(parentDirS,"007_log").getAbsolutePath();

            chrArr = new String[]{"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
            gwin = "0.00001";
            snpWin = "100";
            gridSize = "100000";
            pop1 = "Cultivar";
            pop2 = "Landrace_EU";
            for (int j = 0; j < chrArr.length; j++) {
                int chr = Integer.parseInt(chrArr[j]);
                String pop1fileS = new File(infileDirS,"chr" + chrArr[j] + "_exon_vmap2.1_" + pop1 + "_geno.txt").getAbsolutePath();
                String pop2fileS = new File(infileDirS,"chr" + chrArr[j] + "_exon_vmap2.1_" + pop2 + "_geno.txt").getAbsolutePath();
                String snpInfoS = new File(snpInfoDirS,"chr"+chrArr[j]+"_exon_vmap2.1.snp.txt").getAbsolutePath();
                String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_" + pop1 + "_VS_" + pop2 + "_exonRegion_" + gwin + "_" + snpWin + "_" + gridSize).getAbsolutePath();
                String logS = new File(logDirS,new File(outfileS).getName().split(".gz")[0]+".log.txt").getAbsolutePath();
                System.out.println("XPCLR -xpclr " + pop1fileS + " " + pop2fileS + " " +
                        snpInfoS + " " + outfileS + " -w1 " + gwin + " " + snpWin + " "+
                        gridSize + " " + chr + " -p0 0.95 " + " > " + logS + " &");
            }
        }
        else if (ploidy.equals("ab")){
            infileDirS = "/data4/home/aoyue/vmap2/analysis/030_XPCLR/005_tetraploid/005_genotype";
            outfileDirS = "/data4/home/aoyue/vmap2/analysis/030_XPCLR/005_tetraploid/006_output";
            snpInfoDirS = "/data4/home/aoyue/vmap2/analysis/030_XPCLR/005_tetraploid/002_snp_file";
            logDirS = "/data4/home/aoyue/vmap2/analysis/030_XPCLR/005_tetraploid/007_log";
            chrArr = new String[]{"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};

            gwin = "0.00001";
            snpWin = "100";
            gridSize = "100";
            pop1 = "Domesticated_emmer";
            pop2 = "Wild_emmer";

            for (int j = 0; j < chrArr.length; j++) {
                int chr = Integer.parseInt(chrArr[j]);
                String pop1fileS = new File(infileDirS,"chr" + chrArr[j] + "_exon_vmap2.1_" + pop1 + "_geno.txt").getAbsolutePath();
                String pop2fileS = new File(infileDirS,"chr" + chrArr[j] + "_exon_vmap2.1_" + pop2 + "_geno.txt").getAbsolutePath();
                String snpInfoS = new File(snpInfoDirS,"chr"+chrArr[j]+"_exon_vmap2.1.snp.txt").getAbsolutePath();
                String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_" + pop1 + "_VS_" + pop2 + "_exonRegion_" + gwin + "_" + snpWin + "_" + gridSize).getAbsolutePath();
                String logS = new File(logDirS,new File(outfileS).getName().split(".gz")[0]+".log.txt").getAbsolutePath();
                System.out.println("XPCLR -xpclr " + pop1fileS + " " + pop2fileS + " " +
                        snpInfoS + " " + outfileS + " -w1 " + gwin + " " + snpWin + " "+
                        gridSize + " " + chr + " -p0 0.95" +
                        " > " + logS + " &");
            }
        }

        //-w1 0.0002 200 2000 1 -p1 0.95
        //XPCLR -xpclr /data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/004_genoFile/chr001_Cultivar_geno.txt /data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/004_genoFile/chr001_Landrace_Europe_geno.txt /data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/001_chrposRefAlt/chr001_exon_vmap2.1.pos.Base.txt chr001_CLvsEU_100kbwindow -w1 0.005 600 100000 1 -p0 0.95 > log_chr001_CLvsEU_100kbwindow.txt 2>&1 &
    }

    /**
     * 根据文件夹获取pop的基因型XPCLR格式，多线程运行,适用于六倍体群体
     */
    public void getGenotype_parallele(String infileDirS,String popfileS,String outfileDirS){

        String pop = new File(popfileS).getName().replaceFirst(".txt","");
        List<String> queryTaxal = new AoFile().getStringListwithoutHeader(popfileS,0); //获取pop列表
        Collections.sort(queryTaxal);

        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = new File(outfileDirS,f.getName().split(".vcf")[0]+ "_" + pop + "_geno.txt").getAbsolutePath();
            try{
                BufferedReader br = new AoFile().readFile(infileS);
                BufferedWriter bw = new AoFile().writeFile(outfileS);
                String temp = null;
                List<String> l = new ArrayList<>();
                List<Integer> indexHexa = new ArrayList<>();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("##"))continue;
                    if (temp.startsWith("#C")){
                        l = PStringUtils.fastSplit(temp);
                        for (int i = 9; i < l.size(); i++) {
                            String taxon = l.get(i);
                            int index = Collections.binarySearch(queryTaxal, taxon);
                            if (index > -1) { //当找到列表中的taxa时，写列表中的taxa信息
                                indexHexa.add(i);
                            }
                        }
                        Collections.sort(indexHexa);
//                    System.out.println("Finish find the pop index from vcffile.");
                    }
                    if (!temp.startsWith("#")) {
                        l = PStringUtils.fastSplit(temp);
                        List<String> lGeno = new ArrayList<>();
                        for (int i = 0; i < indexHexa.size(); i++) { //无论有无基因型，都加进去了
                            lGeno.add(l.get(indexHexa.get(i)));
                        }
                        String[] GenoArray = lGeno.toArray(new String[lGeno.size()]);

                        String geno = this.getGenoInfo(GenoArray);
                        bw.write(geno);
                        bw.newLine();
                        cnt++;
                    }
                }
                System.out.println(cnt + " SNP " + new File(infileS).getName() + " is completed at " + outfileS);
                br.close();
                bw.flush();
                bw.close();
            }catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });

    }



    public void mergeTXT(String infileDirS, String outfileDirS){
        String outfileS = new File(outfileDirS,new File(infileDirS).listFiles()[0].getName().split("_exon_")[1]).getAbsolutePath() ;
        AoFile.mergeTxt(infileDirS,outfileS);
    }

    /**
     * 第二步：根据
     */
    public void getSNPdensity_tetraploid(String infileDirS, String outfileDirS){
        int window = 50000;
        int step = 10000;

        //注意没有36号染色体
        String[] chrArr ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};

        for (int j = 0; j < chrArr.length; j++) {
            String infileS = new File(infileDirS,"chr" + chrArr[j] + "_exon_vmap2.1.pos.Base.txt.gz").getAbsolutePath();
            String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_exon_vmap2.1_snpDensity_" + window + "window_" + step + "step.txt.gz").getAbsolutePath();
            this.calDensity(infileS,1,2,window,step,outfileS);
        }

    }

    /**
     * 第二步：根据
     */
    public void getSNPdensity_hexaploid(String infileDirS, String outfileDirS){
//        int window = 50000;
//        int step = 10000;

        int window = 2000;
        int step = 2000;

        //注意没有36号染色体
        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","037","038","039","040","041","042"};

        for (int j = 0; j < chrArr.length; j++) {
            String infileS = new File(infileDirS,"chr" + chrArr[j] + "_exon_vmap2.1.pos.Base.txt.gz").getAbsolutePath();
            String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_exon_vmap2.1_snpDensity_" + window + "window_" + step + "step.txt.gz").getAbsolutePath();
            this.calDensity(infileS,1,2,window,step,outfileS);
        }

    }


    /**
     *
     * 根据science 发表的重组文件，向 snp info 中添加 genetic position
     */
    public void step2_mkSNPfile (String infileDirS, String outfileDirS) {

        String recombinationFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/007_recombination/001_recombination/iwgsc_refseqv1.0_mapping_data_chrID.txt";

        ColumnTable<String> t = new ColumnTable<>(recombinationFileS);
        System.out.println(t.getRowNumber());
        int chrNum = Integer.parseInt(t.getCell(t.getRowNumber()-1, 1)); //获取最后一行第2列的数字，即染色体最大值，这里是42号染色体
        TIntArrayList[] startLists = new TIntArrayList[chrNum]; //42条染色体中，每条染色体的窗口起始位置的集合
        TFloatArrayList[] geneticPosLists = new TFloatArrayList[chrNum]; //42条染色体中，每条染色体的每个窗口对应cross数值集合
        for (int i = 0; i < startLists.length; i++) { //对每个数组内的集合进行初始化
            startLists[i] = new TIntArrayList();
            geneticPosLists[i] = new TFloatArrayList();
        }
        int index = -1;
        for (int i = 0; i < t.getRowNumber(); i++) {
            index = Integer.parseInt(t.getCell(i, 1))-1; //染色体号的索引，即1号染色体索引为0
            startLists[index].add(Integer.parseInt(t.getCell(i, 2))); //将每个Bin的起始位置加入集合中
            geneticPosLists[index].add(Float.parseFloat(t.getCell(i, 3)));
        }
        for (int i = 0; i < startLists.length; i++) {
            startLists[i].sort();
        }

        List<File> fList = IOUtils.getFileListInDirEndsWith(infileDirS, ".txt.gz");
        fList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = new File(outfileDirS,f.getName().replaceFirst("_exon_vmap2.1.pos.Base.txt.gz","_exon_vmap2.1.snp.txt")).getAbsolutePath();
            try {
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String header = br.readLine();
                String temp = null;
                List<String> l = new ArrayList<>();
                int cnt = 0;
                int chrIndex = -1;
                int posIndex = -1;
                int currentPos = -1;
                String geneticPos = null;
                String ID = null;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    cnt++;
                    ID = l.get(0);
                    chrIndex = Integer.parseInt(l.get(1))-1; // 索引1 含有染色体号
                    currentPos = Integer.parseInt(l.get(2)); //索引2 含有位置
                    posIndex = startLists[chrIndex].binarySearch(currentPos); //在 刚刚的起始集合里搜索 index
                    if (posIndex < 0) posIndex = -posIndex-2;
                    if(posIndex == -1){ //说明在第一个数字前面
                        geneticPos = String.format("%.6f", geneticPosLists[chrIndex].get(0)/100);
                        bw.write(ID + "\t" + l.get(1)+ "\t" + geneticPos + "\t" + l.get(2)+ "\t" + l.get(3)+ "\t" + l.get(4));
                        bw.newLine();
                    }
                    else if(posIndex == startLists[chrIndex].size()-1){ //说明在最后一个数字后面
                        geneticPos = String.format("%.6f", geneticPosLists[chrIndex].get(startLists[chrIndex].size()-1)/100);
                        bw.write(ID + "\t" + l.get(1)+ "\t" + geneticPos + "\t" + l.get(2)+ "\t" + l.get(3)+ "\t" + l.get(4));
                        bw.newLine();
                    }
                    else{ //开始进行该点和集合点的比较，找最近值
                        int pos1 = currentPos - startLists[chrIndex].get(posIndex);
                        int pos2 = startLists[chrIndex].get(posIndex+1) - currentPos;
                        if(pos1 > pos2){ //说明离2近
                            geneticPos = String.format("%.6f", geneticPosLists[chrIndex].get(posIndex+1)/100);
                            bw.write(ID + "\t" + l.get(1)+ "\t" + geneticPos + "\t" + l.get(2)+ "\t" + l.get(3)+ "\t" + l.get(4));
                            bw.newLine();
                        }else{ //说明离1近
                            geneticPos = String.format("%.6f", geneticPosLists[chrIndex].get(posIndex)/100);
                            bw.write(ID + "\t" + l.get(1)+ "\t" + geneticPos + "\t" + l.get(2)+ "\t" + l.get(3)+ "\t" + l.get(4));
                            bw.newLine();
                        }

                    }

                }
                br.close();
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    /**
     * 从VCF文件中提取最基本的信息
     */
    public void step1_mkSNPfile(String infileDirS, String outfileDirS){
        this.getIDChrPosRefAlt(infileDirS,outfileDirS);
    }


    public void getIDChrPosRefAlt(String infileDirS, String outfileDirS){
        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        //去除 chr36 染色体,因为没有重组率信息
        for (int i = 0; i < fsList.size(); i++) {
            int chr = Integer.parseInt(fsList.get(i).getName().substring(3,6));
            if (chr == 36) fsList.get(i).delete();
        }
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS, f.getName().split(".vcf")[0] + ".pos.Base.txt.gz").getAbsolutePath();
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String temp = null;
                int cnt = 0;
                bw.write("ID\tChr\tPos\tRef\tAlt\n");
                List<String> l = null;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        continue;
                    }
                    temp = temp.substring(0, 40); //肯定够
                    l = PStringUtils.fastSplit(temp);
                    StringBuilder sb = new StringBuilder();
                    sb = new StringBuilder(l.get(2));
                    sb.append("\t").append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(3)).append("\t").append(l.get(4));
                    bw.write(sb.toString());
                    bw.newLine();
                    if (cnt % 1000000 == 0) {
                        System.out.println("Output " + String.valueOf(cnt) + " SNPs");
                    }
                    cnt++;
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(String.valueOf(cnt) + " SNPs output from " + f.getAbsolutePath());

            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    /**
     * 获取exonVCF的子集
     */
    public void getExonVCFbyPloidy(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/016_exonVCF";
/////////// 六倍体
//        String taxafileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/000_pop_byPloidy/AABBDD.txt";
//        String outfileDirS ="/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/004_hexaploid/000_exonVCF";
//        List<File> fsList = AoFile.getFileListInDir(infileDirS);
//        fsList.parallelStream().forEach(f -> {
//            String infileS = f.getAbsolutePath();
//            String outfileS = new File(outfileDirS,f.getName()).getAbsolutePath();
//            new CalVCF().extractVCF(infileS,outfileS,taxafileS);
//            System.out.println(f.getName() + "\tis completed at " + outfileS);
//        });

/////////// 四倍体
        String taxafileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/000_pop_byPloidy/AABB.txt";
        String outfileDirS ="/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/005_tetraploid/000_exonVCF";
        String[] chrArr ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};

        for (int i = 0; i < chrArr.length; i++) {
            String chr = chrArr[i];
            String infileS = new File(infileDirS,"chr" + chr + "_exon_vmap2.1.vcf.gz").getAbsolutePath();
            String outfileS = new File(outfileDirS, "chr" + chr + "_exon_vmap2.1.vcf.gz").getAbsolutePath();
            new CalVCF().extractVCF(infileS,outfileS,taxafileS);
            System.out.println(chr + "\tis completed at " + outfileS);
        }
    }

    /**
     * 判断 annnotation DB 中的 site 是否都在 exon VCF 中
     * 结果： 全部都在！
     */
    public void checkAnnotationDBisinExonVCF(){

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/018_exonSNPAnnotation";
        String vcffileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/016_exonVCF";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/038_XPCLR/001_ifAnnotationExonVCF";
        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.stream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String exonVCFfileS = new File(vcffileDirS,f.getName().replaceFirst("_SNP_anno.txt.gz","_exon_vmap2.1.vcf.gz")).getAbsolutePath();
                String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + ".txt").getAbsolutePath();

                // 现获取ExonVCF 的pos库
                TIntArrayList posList = CalVCF.extractVCFPos(exonVCFfileS);
                posList.sort();

                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String header = br.readLine();
                String temp = null;
                List<String> l = new ArrayList<>();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    String chr = l.get(1);
                    int pos = Integer.parseInt(l.get(2));
                    int index = posList.binarySearch(pos);
                    if (index < 0){
                        System.out.println(chr + "\t" + pos);
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName() + "\tis completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    //********************************************** 以下 **************
    //********************************************** 之前的程序 *************************************************************************************
    //********************************************** 之前的程序 *************************************************************************************
    //********************************************** 之前的程序 *************************************************************************************
    //********************************************** 之前的程序 *************************************************************************************
    //********************************************** 之前的程序 *************************************************************************************
    //********************************************** 之前的程序 *************************************************************************************

    public void checkTopGeneDistribution(){
        List<String> l = new ArrayList<>();

        try {
//            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/007_GO/001_input/001_GeneID_v2.txt";
            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/002_DEvsWE/002_merge/003_WEvsDE_exonRegion_0.0001_100_50000.xpclr_addHeader_sortbyXPCLR_top0.01_geneList.txt";
            BufferedReader br = new AoFile().readFile(infileS);
            String temp = null;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnt++;
                String chr = temp.substring(7,9);
                l.add(chr);
            }
            br.close();

            AoMath.countCase_fromList_outFile(l);

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 根据得到的topK的XPCLR结果，进行受选择区域的CHR POS的提取，并判断其所在的基因位置
     * 写出2个文件:①chr pos transcript
     * ②基因列表，无表头
     */
    public void getSelectedPos(String infileS){
        String snpAnnoS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/004_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
        //Top K xpclr regions
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/004_merge/001_CLvsEU_exonRegion_0.0001_200_50000_addHeader_sortbyXPCLR_top0.01.xpclr.txt";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/002_DEvsWE/002_merge/003_WEvsDE_exonRegion_0.0001_100_50000.xpclr_addHeader_sortbyXPCLR_top0.01.txt";
                // 受选择区域的位点列表
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/004_merge/001_CLvsEU_exonRegion_0.0001_200_50000_addHeader_sortbyXPCLR_top0.05_transcript.xpclr.txt";
//        String geneList = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/007_GO/001_input/001_GeneID_v2.txt";

        String outfileS = new File(infileS).getAbsolutePath().replaceFirst(".txt","_transcript.txt");
        String geneList = new File(infileS).getAbsolutePath().replaceFirst(".txt","_geneList.txt");

        Set<String> s = new HashSet<>();
        /**
         * 初始化受选择区域的 起始集合 终止集合
         */
        int chrNum = 42;
        int bin = 50000;
        int currentPos = -1;
        int chrIndex = -1;
        RowTable<String> t = new RowTable<>(infileS);
        TIntList[] startLists = new TIntList[chrNum]; //所有的起始位点建立一个集合
        TIntList[] endLists = new TIntList[chrNum]; //所有的终止位点建立一个集合
        for (int i = 0; i < chrNum; i++) { //对list数组进行初始化
            startLists[i] = new TIntArrayList();
            endLists[i] = new TIntArrayList();
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            chrIndex = t.getCellAsInteger(i,0) -1;
            currentPos = t.getCellAsInteger(i,3);
            startLists[chrIndex].add(currentPos);
            endLists[chrIndex].add(currentPos+bin);
        }
        for (int i = 0; i < startLists.length; i++) {
            startLists[i].sort();
            endLists[i].sort();
        }
        System.out.println("Finished building the region list");

        /**
         *  ################################### 写出受选择区域的位点
         */

        new AoFile().readheader(snpAnnoS);

        try{
            t = new RowTable(snpAnnoS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write("CHROM\tPOS\tTranscript");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) {
                int index = t.getCellAsInteger(i, 1) - 1; //染色体号的索引
                int pos = t.getCellAsInteger(i,2);
                String trans = t.getCell(i,10);
                /**
                 * 对该位点进行判断，看是否在选择区域,不在选择区域就忽略不计
                 */
                int posIndex =-1;
                posIndex = startLists[index].binarySearch(pos);
                if (posIndex < 0) {
                    posIndex = -posIndex-2; //确保该位点在起始位点的右边
                }
                if (posIndex < 0) continue; //如果不在起始位点的右边，那么就不在范围内，跳过该位点
                if (pos >= endLists[index].get(posIndex)) continue; //确保在末端位点的前面，若不在，也舍去
                bw.write(index+1 + "\t" + pos + "\t" + trans);
                bw.newLine();
                s.add(trans);
            }
            bw.flush();
            bw.close();

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        try{
            BufferedWriter bw = AoFile.writeFile(geneList);
            String[] gene = s.toArray(new String[s.size()]);
            Arrays.sort(gene);
            for (int i = 0; i < gene.length; i++) {
                bw.write(gene[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        System.out.println("Finished in getting the selected region at " + outfileS + " and " + geneList);
    }




    /**
     * 将XPCLR的结果添加 gene 结果,只添加xpclr输出文件中的pos单个位点的结果
     */
    public void addGeneID_onlyGridPos(){
        String geneHCFileS = "/Users/Aoyue/Documents/Data/wheat/gene/001_geneHC/geneHC.txt";
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/004_merge/001_CLvsEU_exonRegion_0.0001_200_50000_addHeader_sortbyXPCLR_top0.05.xpclr.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/004_merge/001_CLvsEU_exonRegion_0.0001_200_50000_addHeader_sortbyXPCLR_top0.05_addGeneID.xpclr.txt";

        Set<String> transSet = new HashSet<>();
        //先处理gene的表格，建立区间
        Table t = TablesawUtils.readTsv(geneHCFileS);
        System.out.println(t.structure());
        t.sortAscendingOn("Chr","TranStart");
        IntColumn chrColumn = t.intColumn("chr"); //返回一个类 InColumn
        int chrNum = chrColumn.countUnique(); //意思是一共有42条染色体
        TIntList[] startLists = new TIntList[chrNum]; //所有的起始位点建立一个集合
        TIntList[] endLists = new TIntList[chrNum]; //所有的终止位点建立一个集合
        List<String>[] tranLists = new ArrayList[chrNum]; //每条染色体都有一个list
        for (int i = 0; i < chrNum; i++) { //对list数组进行初始化
            startLists[i] = new TIntArrayList();
            endLists[i] = new TIntArrayList();
            tranLists[i] = new ArrayList();
        }
        for (int i = 0; i < t.rowCount(); i++) {
            startLists[Integer.parseInt(t.getString(i, 2))-1].add(Integer.parseInt(t.getString(i, 3)));
            endLists[Integer.parseInt(t.getString(i, 2))-1].add(Integer.parseInt(t.getString(i, 4)));
            tranLists[Integer.parseInt(t.getString(i, 2))-1].add(t.getString(i, 1));
        }


        //再处理要添加列的文件
        try{
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write(br.readLine() + "\tGeneID");
            bw.newLine();
            String temp = null;
//            String temp = br.readLine();

            String trans = null;
            List<String> l = new ArrayList<>();
            int currentPos = -1;
            int posIndex = -1;
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null){
                sb.setLength(0);
                l= PStringUtils.fastSplit(temp);
                int chrIndex = Integer.parseInt(l.get(0));
                currentPos = Integer.parseInt(l.get(3));
                // 对 currentPos所在region内的所有pos进行判断
//                if(currentPos < 50000){ //说明在起始位点
//                    for (int i = 0; i < currentPos; i++) {
//                        int pos = i+1;
//                        posIndex = startLists[chrIndex].binarySearch(pos);
//                        if (posIndex < 0) {
//                            posIndex = -posIndex-2; //确保该位点在起始位点的右边
//                        }
//                        if (posIndex < 0) continue; //如果不在起始位点的右边，那么就不在范围内，跳过该位点
//                        if (pos >= endLists[chrIndex].get(posIndex)) continue; //确保在末端位点的前面，若不在，也舍去
//                        trans = tranLists[chrIndex].get(posIndex);
//                        transSet.add(trans);
//                    }
//
//                }else if (currentPos > 50000){ //说明在中间区域
//                    for (int i = currentPos-50000; i < currentPos; i++) {
//                        int pos = i+1;
//                        posIndex = startLists[chrIndex].binarySearch(pos);
//                        if (posIndex < 0) {
//                            posIndex = -posIndex-2; //确保该位点在起始位点的右边
//                        }
//                        if (posIndex < 0) continue; //如果不在起始位点的右边，那么就不在范围内，跳过该位点
//                        if (pos >= endLists[chrIndex].get(posIndex)) continue; //确保在末端位点的前面，若不在，也舍去
//                        trans = tranLists[chrIndex].get(posIndex);
//                        transSet.add(trans);
//                    }
//                }
                posIndex = startLists[chrIndex].binarySearch(currentPos);
                if (posIndex < 0) {
                    posIndex = -posIndex-2; //确保该位点在起始位点的右边
                }
                if (posIndex < 0) continue; //如果不在起始位点的右边，那么就不在范围内，跳过该位点
                if (currentPos >= endLists[chrIndex].get(posIndex)) continue; //确保在末端位点的前面，若不在，也舍去
                trans = tranLists[chrIndex].get(posIndex);

                sb.append(temp).append("\t").append(trans);
                bw.write(sb.toString());
                bw.newLine();
            }



            List<String> goalTransl = new ArrayList<>(transSet);
            for (int i = 0; i < goalTransl.size(); i++) {
                bw.write(goalTransl.get(i));
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();

            new AoMath().countCaseInGroup(outfileS,7);
        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    /**
     * 将XPCLR的结果添加 gene 结果, 单个pos向后计算50Kb，再进行区域内的位点判断
     */
    public void addGeneID(){
        String geneHCFileS = "/Users/Aoyue/Documents/Data/wheat/gene/001_geneHC/geneHC.txt";
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/004_merge/001_CLvsEU_exonRegion_0.0001_200_50000_addHeader_sortbyXPCLR_top0.05.xpclr.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/007_GO/001_input/GeneID.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/004_merge/001_CLvsEU_exonRegion_0.0001_200_50000_addHeader_sortbyXPCLR_top0.05_addGeneID.xpclr.txt";

        Set<String> transSet = new HashSet<>();
        //先处理gene的表格，建立区间
        Table t = TablesawUtils.readTsv(geneHCFileS);
        System.out.println(t.structure());
        t.sortAscendingOn("Chr","TranStart");
        IntColumn chrColumn = t.intColumn("chr"); //返回一个类 InColumn
        int chrNum = chrColumn.countUnique(); //意思是一共有42条染色体
        TIntList[] startLists = new TIntList[chrNum]; //所有的起始位点建立一个集合
        TIntList[] endLists = new TIntList[chrNum]; //所有的终止位点建立一个集合
        List<String>[] tranLists = new ArrayList[chrNum]; //每条染色体都有一个list
        for (int i = 0; i < chrNum; i++) { //对list数组进行初始化
            startLists[i] = new TIntArrayList();
            endLists[i] = new TIntArrayList();
            tranLists[i] = new ArrayList();
        }
        for (int i = 0; i < t.rowCount(); i++) {
            startLists[Integer.parseInt(t.getString(i, 2))-1].add(Integer.parseInt(t.getString(i, 3)));
            endLists[Integer.parseInt(t.getString(i, 2))-1].add(Integer.parseInt(t.getString(i, 4)));
            tranLists[Integer.parseInt(t.getString(i, 2))-1].add(t.getString(i, 1));
        }


        //再处理要添加列的文件
        try{
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
//            String temp = null;
            String temp = br.readLine();

            String trans = null;
            List<String> l = new ArrayList<>();
            int currentPos = -1;
            int posIndex = -1;
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null){
                sb.setLength(0);
                l=PStringUtils.fastSplit(temp);
                int chrIndex = Integer.parseInt(l.get(0));
                currentPos = Integer.parseInt(l.get(3));
                // 对 currentPos所在region内的所有pos进行判断
                if(currentPos < 50000){ //说明在起始位点
                    for (int i = 0; i < currentPos; i++) {
                        int pos = i+1;
                        posIndex = startLists[chrIndex].binarySearch(pos);
                        if (posIndex < 0) {
                            posIndex = -posIndex-2; //确保该位点在起始位点的右边
                        }
                        if (posIndex < 0) continue; //如果不在起始位点的右边，那么就不在范围内，跳过该位点
                        if (pos >= endLists[chrIndex].get(posIndex)) continue; //确保在末端位点的前面，若不在，也舍去
                        trans = tranLists[chrIndex].get(posIndex);
                        transSet.add(trans);
                    }

                }else if (currentPos > 50000){ //说明在中间区域
                    for (int i = currentPos-50000; i < currentPos; i++) {
                        int pos = i+1;
                        posIndex = startLists[chrIndex].binarySearch(pos);
                        if (posIndex < 0) {
                            posIndex = -posIndex-2; //确保该位点在起始位点的右边
                        }
                        if (posIndex < 0) continue; //如果不在起始位点的右边，那么就不在范围内，跳过该位点
                        if (pos >= endLists[chrIndex].get(posIndex)) continue; //确保在末端位点的前面，若不在，也舍去
                        trans = tranLists[chrIndex].get(posIndex);
                        transSet.add(trans);
                    }
                }
            }

            List<String> goalTransl = new ArrayList<>(transSet);
            for (int i = 0; i < goalTransl.size(); i++) {
                bw.write(goalTransl.get(i));
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();

//            new AoMath().countCaseInGroup(outfileS,7);
        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    //获取topK的结果，只输出CHROM POS	Genetic_pos	XPCLR_score
    public void getTopK(String infileS){
        double k = 0.01;
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/004_merge/001_CLvsEU_exonRegion_0.0001_200_50000_addHeader_sortbyXPCLR.xpclr.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/004_merge/001_CLvsEU_exonRegion_0.0001_200_50000_addHeader_sortbyXPCLR_top0.05.xpclr.txt";
        String outfileS = new File(infileS).getAbsolutePath().replaceFirst(".txt","_top" + k + ".txt");
        try{
            int n = AoFile.countFileRowNumber(infileS);
            double line = k*n;
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write(br.readLine());
            bw.newLine();
            for (int i = 0; i < line; i++) {
                bw.write(br.readLine());
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println("Finished in getting top" + k + " at " + outfileS);
            System.out.println("-----------------------------------------------");
            System.out.println("Begin to get the selected region pos");
            this.getSelectedPos(outfileS);
        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    public void sortbyXPCLR(String infileS){
        double xpclr = Double.MIN_VALUE;
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/004_merge/001_CLvsEU_exonRegion_0.0001_200_50000_addHeader.xpclr.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/004_merge/001_CLvsEU_exonRegion_0.0001_200_50000_addHeader_sortbyXPCLR.xpclr.txt";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/002_DEvsWE/002_merge/003_WEvsDE_exonRegion_0.0001_100_50000.xpclr_addHeader.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/002_DEvsWE/002_merge/003_WEvsDE_exonRegion_0.0001_100_50000.xpclr_addHeader_sortbyXPCLR.txt";

        try{
            BufferedReader br = AoFile.readFile(infileS);
            String outfileS = new File(infileS).getAbsolutePath().replaceFirst(".txt","_sortbyXPCLR.txt");
            BufferedWriter bw = AoFile.writeFile(outfileS);

            //进行 list 的构建
            RowTable<String> t = new RowTable<>(infileS);
            List<Record> l = new ArrayList<>();
            for (int i = 0; i < t.getRowNumber(); i++) {
                String chr = t.getCell(i,0);
                int grid = t.getCellAsInteger(i,1);
                int N_SNPs = t.getCellAsInteger(i,2);
                int pos = t.getCellAsInteger(i,3);
                double genetic_pos = t.getCellAsDouble(i,4);
                String xpclrS = t.getCell(i,5);
                if (xpclrS.equals("inf")){
                    xpclr=0.0;
                }else {
                    xpclr = Double.parseDouble(xpclrS);
                }
                if (t.getCell(i,6).isEmpty())continue;
                double max_s = t.getCellAsDouble(i,6);
                Record o = new Record(chr,grid,N_SNPs,pos,genetic_pos,xpclr,max_s);
                l.add(o);
            }

            //排序并写出排序后的文件
            Collections.sort(l);
            bw.write(br.readLine());
            bw.newLine();
            for (int i = 0; i < l.size(); i++) {
                bw.write(l.get(i).getString());
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println("Finished in sorting the XP-CLR result at " + outfileS);
            System.out.println("-----------------------------------------------------");
            System.out.println("Begin to get TopK result");
            this.getTopK(outfileS);

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }


    class Record implements Comparable<Record>{
        public String chr;
        public int grid;
        public int N_SNPs;
        public int pos;
        public double genetic_pos;
        public double xpclr;
        public double max_s;

        public Record(String chr, int grid, int N_SNPs, int pos, double genetic_pos, double xpclr,double max_s){
            this.chr=chr;
            this.grid=grid;
            this.N_SNPs=N_SNPs;
            this.pos=pos;
            this.genetic_pos=genetic_pos;
            this.xpclr=xpclr;
            this.max_s=max_s;
        }

        public String getString(){
            StringBuilder sb = new StringBuilder();
            sb.append(this.chr).append("\t").append(this.grid).append("\t").append(this.N_SNPs).append("\t").
                    append(this.pos).append("\t").append(genetic_pos).append("\t").append(this.xpclr).
                    append("\t").append(this.max_s);
            return sb.toString();
        }

        @Override
        public int compareTo(Record o){
            if(this.xpclr < o.xpclr){
                return 1; //从高到低排序
            }else if (this.xpclr == o.xpclr){
                return 0;
            }
            else{
                return -1;
            }
        }
    }

    /**
     * 将结果不进行坐标转换，只添加表头，把信息不完全的行删除
     *
     */
    public void convertXPCLRCoordinate2(){

        HashMap<String,Integer> hm = new HashMap<String, Integer>();
        String[] chrArr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
        for (int i = 0; i < chrArr.length; i++) {
            String chr = chrArr[i];
            hm.put(chr,i+1);
        }
        try {
//            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/002_merge/CLvsEU_exonRegion_100kbwindow.xpclr.txt";
//            String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/002_merge/001_CLvsEU_exonRegion_100kbwindow_changeChrPos.xpclr.txt";

//            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/004_merge/CLvsEU_exonRegion_0.0001_200_50000.xpclr.txt";
//            String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/004_merge/001_CLvsEU_exonRegion_0.0001_200_50000_addHeader.xpclr.txt";

            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/002_DEvsWE/002_merge/001_WEvsDE_exonRegion_0.0001_100_50000.xpclr.txt";
            String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/002_DEvsWE/002_merge/003_WEvsDE_exonRegion_0.0001_100_50000.xpclr_addHeader.txt";
            BufferedReader br = new AoFile().readFile(infileS);
            BufferedWriter bw = new AoFile().writeFile(outfileS);
            bw.write("CHROM\tGrid\tN_SNPs\tPOS\tGenetic_pos\tXPCLR_score\tMax_s");
            bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp," ");
                if (l.size() < 7)continue;
                cnt++;
                //2 217 1152 21716707.000000 0.651300 584.597428 0.000000
                String chrS = l.get(0);
                String posS = l.get(3);
                int chrID = Integer.parseInt(chrS);
                int posID = (int) Double.parseDouble(posS);
//                String Chr = RefV1Utils.getChromosome(chrID,posID);
//                int pos = RefV1Utils.getPosOnChromosome(chrID,posID);
//                int ID = hm.get(Chr);
                StringBuilder sb = new StringBuilder();
                sb.append(chrID).append("\t").append(l.get(1)).append("\t").append(l.get(2)).append("\t").
                        append(posID).append("\t").append(l.get(4)).append("\t").append(l.get(5)).
                        append("\t").append(l.get(6));
                bw.write(sb.toString());
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println("Finished in mergeing the origin XPCLR output and add header at " + outfileS);
            System.out.println("------------------------------------------------------------");
            System.out.println("Begin to sort XP-CLR result.");
            this.sortbyXPCLR(outfileS);

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 将结果进行坐标转换，并添加表头
     *
     */
    public void convertXPCLRCoordinate(){

        HashMap<String,Integer> hm = new HashMap<String, Integer>();
        String[] chrArr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
        for (int i = 0; i < chrArr.length; i++) {
            String chr = chrArr[i];
            hm.put(chr,i+1);
        }
        try {
//            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/002_merge/CLvsEU_exonRegion_100kbwindow.xpclr.txt";
//            String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/002_merge/001_CLvsEU_exonRegion_100kbwindow_changeChrPos.xpclr.txt";

//            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/002_merge/CLvsEU_exonRegion_0.005m_1400snp_100kbwindoww.xpclr.txt";
//            String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/002_merge/002_CLvsEU_exonRegion_0.005m_1400snp_100kbwindoww_changeChrPos.xpclr.txt";

//            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/004_merge/CLvsEU_exonRegion_0.0001_200_50000.xpclr.txt";
//            String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/004_merge/002_CLvsEU_exonRegion_0.0001_200_50000_changeChrPos_addChr007.xpclr.txt";

            String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/002_DEvsWE/002_merge/001_WEvsDE_exonRegion_0.0001_100_50000.xpclr.txt";
            String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/002_DEvsWE/002_merge/002_WEvsDE_exonRegion_0.0001_100_50000.xpclr_changeChrPos.txt";

            BufferedReader br = new AoFile().readFile(infileS);
            BufferedWriter bw = new AoFile().writeFile(outfileS);
            bw.write("CHROM\tGrid\tN_SNPs\tPOS\tGenetic_pos\tXPCLR_score\tMax_s\tID_Manhatton");
            bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp," ");
                if (l.size() < 7)continue; //如果该行的list小于7，说明没有读写完全就终止了，是最后一行，应该去除
                cnt++;
                //2 217 1152 21716707.000000 0.651300 584.597428 0.000000
                String chrS = l.get(0);
                String posS = l.get(3);
                int chrID = Integer.parseInt(chrS);
                int posID = (int) Double.parseDouble(posS);
                String Chr = RefV1Utils.getChromosome(chrID,posID);
                int pos = RefV1Utils.getPosOnChromosome(chrID,posID);
                int ID = hm.get(Chr);
                StringBuilder sb = new StringBuilder();
                sb.append(Chr).append("\t").append(l.get(1)).append("\t").append(l.get(2)).append("\t").
                        append(pos).append("\t").append(l.get(4)).append("\t").append(l.get(5)).
                        append("\t").append(l.get(6)).append("\t").append(ID);
                bw.write(sb.toString());
                bw.newLine();
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

    public void mergeTxt2(){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/001_out";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/002_merge/CLvsEU_exonRegion_0.005m_1400snp_100kbwindoww.xpclr.txt";
//        new AoFile().mergeTxtwithoutHeader(infileDirS,outfileS);

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/003_out_snpWin200";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/004_merge/CLvsEU_exonRegion_0.0001_200_50000.xpclr.txt";
//        new AoFile().mergeTxtwithoutHeader(infileDirS,outfileS);

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/002_DEvsWE/001_out";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/002_DEvsWE/002_merge/WEvsDE_exonRegion_0.0001_100_50000.xpclr.txt";
        AoFile.mergeTxtwithoutHeader(infileDirS,outfileS);

    }

    /**
     *
     * 以 100kb为窗口，统计各个亚基因组的SNP个数是多少
     */
    public void statisticSNPdensity(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/004_snpDensity/exon_vmap2.1.pos.Base.density_100k.txt";
        String hmfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_chrList/ChrID.txt";
        HashMap<String,String> hm = new AoFile().getHashMapStringKey(hmfileS,0,3);
        new AoFile().addColumbyString(infileS,0,hm,"Subgenome");
    }

    public void mergeTxt(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/004_snpDensity/001_calDensity_100kb";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/004_snpDensity/exon_vmap2.1.pos.Base.density_100k.txt";
        new CountSites().mergeTxt(infileDirS,outfileS);
    }

    public void script_calSNPdensity(){
        String infileDirS = "/data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/001_chrposRefAlt";
        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/003_calDensity";
        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
        for (int j = 0; j < chrArr.length; j++) {
            String infileS = new File(infileDirS,"chr" + chrArr[j] + "_exon_vmap2.1.pos.Base.txt").getAbsolutePath();
            String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_exon_vmap2.1.pos.Base.density_100k.txt").getAbsolutePath();
            System.out.println("java -jar 039_calDensity.jar " + infileS + " 1 3 100000 100000 " + outfileS + " &");
        }

    }

    public void script_XPCLR_tetraploid(){
//        XPCLR -xpclr /data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/004_genoFile/chr001_Cultivar_geno.txt /data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/004_genoFile/chr001_Landrace_Europe_geno.txt /data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/001_chrposRefAlt/chr001_exon_vmap2.1.pos.Base.txt chr001_CLvsEU_100kbwindow -w1 0.005 600 100000 1 -p0 0.95 > log_chr001_CLvsEU_100kbwindow.txt 2>&1 &
        String infileDirS = "/data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/004_genoFile";
//        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/008_out_tetraploid_snpWin100";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/009_out_tetraploid_snpWin100";
        String snpInfoDirS = "/data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/001_chrposRefAlt";
        String logDirS = "/data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/log";
        String[] chrArr ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};

//        //在203服务器上运行试试
//
//        String infileDirS = "/data1/home/aoyue/vmap2/analysis/001_XPCLR/001_exon/002_genoFile";
//        String outfileDirS = "/data1/home/aoyue/vmap2/analysis/001_XPCLR/001_exon/003_out_tetraploid_snpWin100";
//        String snpInfoDirS = "/data1/home/aoyue/vmap2/analysis/001_XPCLR/001_exon/001_snpInfo";
//        String logDirS = "/data1/home/aoyue/vmap2/analysis/001_XPCLR/001_exon/log";
//        String[] chrArr ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};


        String gwin = "0.0001";
        String snpWin = "100";
        String gridSize = "50000";

        for (int j = 0; j < chrArr.length; j++) {
            int chr = Integer.parseInt(chrArr[j]);
            String snpInfoS = new File(snpInfoDirS,"chr"+chrArr[j]+"_exon_vmap2.1.pos.Base.txt").getAbsolutePath();

//            String pop1fileS = new File(infileDirS,"chr" + chrArr[j] + "_exon_vmap2.1_Domesticated_emmer_geno.txt").getAbsolutePath();
//            String pop2fileS = new File(infileDirS,"chr" + chrArr[j] + "_exon_vmap2.1_Wild_emmer_geno.txt").getAbsolutePath();
//            String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_DEvsFTT_exonRegion_"+gwin+"_"+snpWin+"_"+gridSize).getAbsolutePath();

            String pop1fileS = new File(infileDirS,"chr" + chrArr[j] + "_exon_vmap2.1_Free_threshing_tetraploid_geno.txt").getAbsolutePath();
            String pop2fileS = new File(infileDirS,"chr" + chrArr[j] + "_exon_vmap2.1_Domesticated_emmer_geno.txt").getAbsolutePath();
            String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_FTTvsDE_exonRegion_"+gwin+"_"+snpWin+"_"+gridSize).getAbsolutePath();

            String logS = new File(logDirS,new File(outfileS).getName().split(".gz")[0]+".txt").getAbsolutePath();
            System.out.println("XPCLR -xpclr " + pop1fileS + " " + pop2fileS + " " +
                    snpInfoS + " " + outfileS + " -w1 " + gwin + " " + snpWin + " "+
                    gridSize + " " + chr + " -p0 0.95" +
                    " > " + logS);
        }

    }


    public void script_XPCLR(){
        //-w1 0.0002 200 2000 1 -p1 0.95
        //XPCLR -xpclr /data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/004_genoFile/chr001_Cultivar_geno.txt /data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/004_genoFile/chr001_Landrace_Europe_geno.txt /data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/001_chrposRefAlt/chr001_exon_vmap2.1.pos.Base.txt chr001_CLvsEU_100kbwindow -w1 0.005 600 100000 1 -p0 0.95 > log_chr001_CLvsEU_100kbwindow.txt 2>&1 &
        String infileDirS = "/data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/004_genoFile";
//        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/005_out";
        String outfileDirS = "/data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/007_out_snpWin200";
        String snpInfoDirS = "/data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/001_chrposRefAlt";
        String logDirS = "/data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/log";
        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
//        String gwin = "0.005";
//        String gwin = "0.0002";
        String gwin = "0.0001";
//        String snpWin = "1400";
        String snpWin = "200";
//        String snpWin = "100";
        String gridSize = "50000";
//        String gridSize = "100000";
//        String gridSize = "2000";

        for (int j = 0; j < chrArr.length; j++) {
            int chr = Integer.parseInt(chrArr[j]);
            String pop1fileS = new File(infileDirS,"chr" + chrArr[j] + "_exon_vmap2.1_Cultivar_geno.txt").getAbsolutePath();
            String pop2fileS = new File(infileDirS,"chr" + chrArr[j] + "_exon_vmap2.1_Landrace_Europe_geno.txt").getAbsolutePath();
            String snpInfoS = new File(snpInfoDirS,"chr"+chrArr[j]+"_exon_vmap2.1.pos.Base.txt").getAbsolutePath();
            String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_CLvsEU_exonRegion_"+gwin+"_"+snpWin+"_"+gridSize).getAbsolutePath();
            String logS = new File(logDirS,new File(outfileS).getName().split(".gz")[0]+".txt").getAbsolutePath();
            System.out.println("XPCLR -xpclr " + pop1fileS + " " + pop2fileS + " " +
                    snpInfoS + " " + outfileS + " -w1 " + gwin + " " + snpWin + " "+
                    gridSize + " " + chr + " -p0 0.95" +
                    " > " + logS);
        }

    }

    /**
     * 根据文件夹获取pop的基因型XPCLR格式，多线程运行,适用于四倍体群体
     */
    public void getGenotypeXPCLR_parallele_diploid(String infileDirS,String popfileS,String outfileDirS){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/022_subsetVCF/001_singleChr0.001";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/003_genoTest/002_geno";
//        String popfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/003_genoTest/001_pop/Cultivar.txt";
        String pop = new File(popfileS).getName().replaceFirst(".txt","");
        List<String> queryTaxal = new AoFile().getStringListwithoutHeader(popfileS,0); //获取pop列表
        String[] chrArr ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};
        Arrays.sort(chrArr);
        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String chr = f.getName().substring(3,6);
            int jj = Arrays.binarySearch(chrArr,chr);
            if (jj > -1){
                String outfileS = new File(outfileDirS,f.getName().split(".vcf")[0]+ "_" + pop + "_geno.txt").getAbsolutePath();
                try{
                    BufferedReader br = new AoFile().readFile(infileS);
                    BufferedWriter bw = new AoFile().writeFile(outfileS);
                    String temp = null;
                    List<String> l = new ArrayList<>();
                    List<Integer> indexHexa = new ArrayList<>();
                    int cnt = 0;
                    while ((temp = br.readLine()) != null) {
                        if (temp.startsWith("##"))continue;
                        if (temp.startsWith("#C")){
                            l = PStringUtils.fastSplit(temp);
                            for (int i = 9; i < l.size(); i++) {
                                String taxon = l.get(i);
                                int index = Collections.binarySearch(queryTaxal, taxon);
                                if (index > -1) { //当找到列表中的taxa时，写列表中的taxa信息
                                    indexHexa.add(i);
                                }
                            }
                            Collections.sort(indexHexa);
//                    System.out.println("Finish find the pop index from vcffile.");
                        }
                        if (!temp.startsWith("#")) {
                            l = PStringUtils.fastSplit(temp);
                            List<String> lGeno = new ArrayList<>();
                            for (int i = 0; i < indexHexa.size(); i++) { //无论有无基因型，都加进去了
                                lGeno.add(l.get(indexHexa.get(i)));
                            }
                            String[] GenoArray = lGeno.toArray(new String[lGeno.size()]);

                            String geno = this.getGenoInfo(GenoArray);
                            bw.write(geno);
                            bw.newLine();
                            cnt++;
                        }
                    }
                    System.out.println(cnt + " SNP " + new File(infileS).getName() + " is completed at " + outfileS);
                    br.close();
                    bw.flush();
                    bw.close();
                }catch (Exception e) {
                    e.printStackTrace();
                    System.exit(1);
                }

            }

        });

    }

    /**
     * 根据文件夹获取pop的基因型XPCLR格式，多线程运行,适用于四倍体群体
     */
    public void getGenotypeXPCLR_parallele_tetra(String infileDirS,String popfileS,String outfileDirS){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/022_subsetVCF/001_singleChr0.001";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/003_genoTest/002_geno";
//        String popfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/003_genoTest/001_pop/Cultivar.txt";
        String pop = new File(popfileS).getName().replaceFirst(".txt","");
        List<String> queryTaxal = new AoFile().getStringListwithoutHeader(popfileS,0); //获取pop列表
        String[] chrArr ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};
        Arrays.sort(chrArr);
        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String chr = f.getName().substring(3,6);
            int jj = Arrays.binarySearch(chrArr,chr);
            if (jj > -1){
                String outfileS = new File(outfileDirS,f.getName().split(".vcf")[0]+ "_" + pop + "_geno.txt").getAbsolutePath();
                try{
                    BufferedReader br = new AoFile().readFile(infileS);
                    BufferedWriter bw = new AoFile().writeFile(outfileS);
                    String temp = null;
                    List<String> l = new ArrayList<>();
                    List<Integer> indexHexa = new ArrayList<>();
                    int cnt = 0;
                    while ((temp = br.readLine()) != null) {
                        if (temp.startsWith("##"))continue;
                        if (temp.startsWith("#C")){
                            l = PStringUtils.fastSplit(temp);
                            for (int i = 9; i < l.size(); i++) {
                                String taxon = l.get(i);
                                int index = Collections.binarySearch(queryTaxal, taxon);
                                if (index > -1) { //当找到列表中的taxa时，写列表中的taxa信息
                                    indexHexa.add(i);
                                }
                            }
                            Collections.sort(indexHexa);
//                    System.out.println("Finish find the pop index from vcffile.");
                        }
                        if (!temp.startsWith("#")) {
                            l = PStringUtils.fastSplit(temp);
                            List<String> lGeno = new ArrayList<>();
                            for (int i = 0; i < indexHexa.size(); i++) { //无论有无基因型，都加进去了
                                lGeno.add(l.get(indexHexa.get(i)));
                            }
                            String[] GenoArray = lGeno.toArray(new String[lGeno.size()]);

                            String geno = this.getGenoInfo(GenoArray);
                            bw.write(geno);
                            bw.newLine();
                            cnt++;
                        }
                    }
                    System.out.println(cnt + " SNP " + new File(infileS).getName() + " is completed at " + outfileS);
                    br.close();
                    bw.flush();
                    bw.close();
                }catch (Exception e) {
                    e.printStackTrace();
                    System.exit(1);
                }

            }

        });

    }

    /**
     * 根据文件夹获取pop的基因型XPCLR格式，多线程运行,适用于六倍体群体
     */
    public void getGenotypeXPCLR_parallele(String infileDirS,String popfileS,String outfileDirS){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/022_subsetVCF/001_singleChr0.001";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/003_genoTest/002_geno";
//        String popfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/003_genoTest/001_pop/Cultivar.txt";
        String pop = new File(popfileS).getName().replaceFirst(".txt","");
        List<String> queryTaxal = new AoFile().getStringListwithoutHeader(popfileS,0); //获取pop列表

        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = new File(outfileDirS,f.getName().split(".vcf")[0]+ "_" + pop + "_geno.txt").getAbsolutePath();
            try{
                BufferedReader br = new AoFile().readFile(infileS);
                BufferedWriter bw = new AoFile().writeFile(outfileS);
                String temp = null;
                List<String> l = new ArrayList<>();
                List<Integer> indexHexa = new ArrayList<>();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("##"))continue;
                    if (temp.startsWith("#C")){
                        l = PStringUtils.fastSplit(temp);
                        for (int i = 9; i < l.size(); i++) {
                            String taxon = l.get(i);
                            int index = Collections.binarySearch(queryTaxal, taxon);
                            if (index > -1) { //当找到列表中的taxa时，写列表中的taxa信息
                                indexHexa.add(i);
                            }
                        }
                        Collections.sort(indexHexa);
//                    System.out.println("Finish find the pop index from vcffile.");
                    }
                    if (!temp.startsWith("#")) {
                        l = PStringUtils.fastSplit(temp);
                        List<String> lGeno = new ArrayList<>();
                        for (int i = 0; i < indexHexa.size(); i++) { //无论有无基因型，都加进去了
                            lGeno.add(l.get(indexHexa.get(i)));
                        }
                        String[] GenoArray = lGeno.toArray(new String[lGeno.size()]);

                        String geno = this.getGenoInfo(GenoArray);
                        bw.write(geno);
                        bw.newLine();
                        cnt++;
                    }
                }
                System.out.println(cnt + " SNP " + new File(infileS).getName() + " is completed at " + outfileS);
                br.close();
                bw.flush();
                bw.close();
            }catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });

    }

    /**
     *
     * java运行，不生成log文件
     */
    public void script_getGenotypeXPCLR2(){
//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/011_VMapII";
//        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/022_XPCLR/001_prepare/004_genoFile";

        String infileDirS = "/data4/home/aoyue/vmap2/feilu/002_genicSNP/002_exonSNPVCF";
        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/004_genoFile";

        String popfileS = "/data4/home/aoyue/vmap2/analysis/022_XPCLR/001_prepare/000_pop/Cultivar.txt";
//        String popfileS = "/data4/home/aoyue/vmap2/analysis/022_XPCLR/001_prepare/000_pop/Landrace_Europe.txt";

        String pop = new File(popfileS).getName().replaceFirst(".txt","");
        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
        for (int j = 0; j < chrArr.length; j++) {
            String infileS = new File(infileDirS,"chr" + chrArr[j] + "_exon_vmap2.1.vcf.gz").getAbsolutePath();
            String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_"+pop+"_geno.txt").getAbsolutePath();
            System.out.println("java -jar 040_getGenotypeXPCLR.jar " + infileS + " " + popfileS + " " + outfileS );
        }

    }

    /**
     *
     * java运行，一个java生成一个log文件
     */
    public void script_getGenotypeXPCLR(){
//        String infileDirS = "/data4/home/aoyue/vmap2/genotype/mergedVCF/011_VMapII";
//        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/022_XPCLR/001_prepare/004_genoFile";

        String infileDirS = "/data4/home/aoyue/vmap2/feilu/002_genicSNP/002_exonSNPVCF";
        String outfileDirS ="/data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/004_genoFile";

        String logDirS = "/data4/home/aoyue/vmap2/aaPlantGenetics/log_040";
        String popfileS = "/data4/home/aoyue/vmap2/analysis/022_XPCLR/001_prepare/000_pop/Cultivar.txt";
//        String popfileS = "/data4/home/aoyue/vmap2/analysis/022_XPCLR/001_prepare/000_pop/Landrace_Europe.txt";

        String pop = new File(popfileS).getName().replaceFirst(".txt","");
        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
        for (int j = 0; j < chrArr.length; j++) {
            String infileS = new File(infileDirS,"chr" + chrArr[j] + "_exon_vmap2.1.vcf.gz").getAbsolutePath();
            String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_"+pop+"_geno.txt").getAbsolutePath();
            String logfileS = new File(logDirS,"log_" + new File(outfileS).getName().split(".gz")[0]).getAbsolutePath(); //不管是不是gz结尾，我们只取gz前的部分，妙！
            System.out.println("nohup java -jar 040_getGenotypeXPCLR.jar " + infileS + " " + popfileS + " " + outfileS + " > " + logfileS  + " 2>&1 &" );
        }

    }

    /**
     *
     *
     */
    public void getGenotypeXPCLR(String infileS, String popfileS, String outfileS){
//        String popfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/003_genoTest/001_pop/Cultivar.txt";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/022_subsetVCF/001_singleChr0.001/chr001.subgenome.maf0.01.SNP_bi.subset.vcf.gz";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/003_genoTest/002_geno/chr001_cultivar.geno.txt";
        List<String> queryTaxal = new AoFile().getStringListwithoutHeader(popfileS,0); //获取pop列表

        try {
            BufferedReader br = new AoFile().readFile(infileS);
            BufferedWriter bw = new AoFile().writeFile(outfileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            List<Integer> indexHexa = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##"))continue;
                if (temp.startsWith("#C")){
                    l = PStringUtils.fastSplit(temp);
                    for (int i = 9; i < l.size(); i++) {
                        String taxon = l.get(i);
                        int index = Collections.binarySearch(queryTaxal, taxon);
                        if (index > -1) { //当找到列表中的taxa时，写列表中的taxa信息
                            indexHexa.add(i);
                        }
                    }
                    Collections.sort(indexHexa);
//                    System.out.println("Finish find the pop index from vcffile.");
                }
                if (!temp.startsWith("#")) {
                    l = PStringUtils.fastSplit(temp);
                    List<String> lGeno = new ArrayList<>();
                    for (int i = 0; i < indexHexa.size(); i++) { //无论有无基因型，都加进去了
                        lGeno.add(l.get(indexHexa.get(i)));
                    }
                    String[] GenoArray = lGeno.toArray(new String[lGeno.size()]);

                    String geno = this.getGenoInfo(GenoArray);
                    bw.write(geno);
                    bw.newLine();
                    cnt++;
                }
            }
            System.out.println(cnt + " SNP " + new File(infileS).getName() + " is completed at " + outfileS);
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    public String getGenoInfo(String[] genoArray) {
        String geno = null;
        List<String> tempList = null;
        List<String> temList = null;
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < genoArray.length; i++) {
            tempList = PStringUtils.fastSplit(genoArray[i], ":"); //tempList是包含基因型GT AD PL的集合
            temList = PStringUtils.fastSplit(tempList.get(0), "/"); //temList是包含基因型拆分后的集合
            String ref1 = temList.get(0);
            String alt2 = temList.get(1);
            if (ref1.equals(".")){
                ref1 = "9";
                alt2 = "9";
            }
            sb.append(ref1).append(" ").append(alt2).append(" ");
        }
        geno = sb.toString();
        return geno;
    }


    /**
     * 根据 chr pos 两列,返回每个window内的变异个数
     * @param infileS
     * @param chrIndex 染色体所在的那一列
     * @param posIndex 位置所在的那一列
     * @param window 滑窗大小
     * @param step 步移大小
     * @param outfileS
     */
    public void calDensity(String infileS,int chrIndex,int posIndex, int window, int step, String outfileS){
        // 1.将pos转为为list,找到最大值，根据最大值确定bin的数目
        // 2.建立count数目和list数组，对pos进行循环，找到每个bin的左边的数目和value的集合，求这个集合的平均值，最大值，方差等
        // 3.输出，每个bin的值
        //************************

        //先求最大值,即染色体长度
        RowTable<String> t = new RowTable(infileS);
        int chrlength = Integer.valueOf(t.getCell(t.getRowNumber() - 1, posIndex)); //t.getRowNumber()是文件的行数，不包括header。 这里getCell得到的是索引         //染色体的长度是最后一行pos的位置
        System.out.println(new File(infileS).getName().substring(0,5) + " length is " + chrlength);
        String chr = t.getCell(0,chrIndex);

        int[][] bound = this.initializeWindowStep(chrlength, window,step);
        int[] count = new int[bound.length]; //查看每个bin里面的变异个数
        int[] boundright = new int[bound.length]; //只看左边的bound
        int[] boundleft = new int[bound.length]; //右边的bound
        for (int i = 0; i < bound.length; i++) { //每个bound的左边
            boundleft[i] = bound[i][0];
            boundright[i] = bound[i][1];
        }

        for (int i = 0; i < t.getRowNumber(); i++) {
            int pos = Integer.parseInt(t.getCell(i,posIndex));

            int indexleft = Arrays.binarySearch(boundleft, pos);
            if (indexleft < 0) {
                indexleft = -indexleft - 2 +1;
            }
            else if (indexleft > -1) {
                indexleft = indexleft +1;
            }


            int indexright = Arrays.binarySearch(boundright, pos);
            if (indexright < 0){
                indexright = -indexright-1;
            }
            else if (indexright > -1){
                indexright = indexright +1;
            }


            for (int j = indexright; j < indexleft ; j++) {
                count[j]++; //每个Bin 里面的变异个数
            }
        }

        try {
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write("CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tSub");
            bw.newLine();
            for (int i = 0; i < bound.length; i++) {
                StringBuilder sb = new StringBuilder();
                String sub = AoString.getSubfromChrID(Integer.parseInt(chr));
                sb.append(chr).append("\t").append(bound[i][0]).append("\t").append(bound[i][1]-1).append("\t").append(count[i]).append("\t").append(sub);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println( "Bin calculation is completed at "+ outfileS);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }


    /**
     * return the bound needed
     *
     * @param chrLength
     * @param windowSize
     * @param windowStep
     * @return
     */
    private int[][] initializeWindowStep (int chrLength, int windowSize, int windowStep) {

        TIntArrayList startList = new TIntArrayList();
        TIntArrayList endList = new TIntArrayList();
        int start = 1;
        int end = start+windowSize;
        while (end < chrLength) { //这里将 start < chrLength 改为 end < chrLength
            startList.add(start);
            endList.add(end);
            start+=windowStep;
            end = start+windowSize;
        }

        int[][] bound = new int[startList.size()][2];
        for (int i = 0; i < startList.size(); i++) {
            bound[i][0] = startList.get(i);
            bound[i][1] = endList.get(i);
        }
        return bound;
    }


    /**
     *
     * 根据science 发表的重组文件，向 snp info 中添加 genetic position
     */
    public void addgeneticPos () {
//        String dirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/002_snp";
//        String recombinationFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/007_recombination/001_recombination/iwgsc_refseqv1.0_mapping_data_chrID.txt";

//        String dirS = "/data4/home/aoyue/vmap2/analysis/022_XPCLR/001_prepare/001_chrposRefAlt";
//        String recombinationFileS = "/data4/home/aoyue/vmap2/analysis/022_XPCLR/001_prepare/iwgsc_refseqv1.0_mapping_data_chrID.txt";

        String dirS = "/data4/home/aoyue/vmap2/analysis/022_XPCLR/002_exon/001_chrposRefAlt";
        String recombinationFileS = "/data4/home/aoyue/vmap2/analysis/022_XPCLR/001_prepare/iwgsc_refseqv1.0_mapping_data_chrID.txt";

        ColumnTable<String> t = new ColumnTable<>(recombinationFileS);
        int chrNum = Integer.parseInt(t.getCell(t.getRowNumber()-1, 1)); //获取最后一行第0列的数字，即染色体最大值，这里是42号染色体
        TIntArrayList[] startLists = new TIntArrayList[chrNum]; //42条染色体中，每条染色体的窗口起始位置的集合
        TFloatArrayList[] geneticPosLists = new TFloatArrayList[chrNum]; //42条染色体中，每条染色体的每个窗口对应cross数值集合
        for (int i = 0; i < startLists.length; i++) { //对每个数组内的集合进行初始化
            startLists[i] = new TIntArrayList();
            geneticPosLists[i] = new TFloatArrayList();
        }
        int index = -1;
        for (int i = 0; i < t.getRowNumber(); i++) {
            index = Integer.parseInt(t.getCell(i, 1))-1; //染色体号的索引，即1号染色体索引为0
            startLists[index].add(Integer.parseInt(t.getCell(i, 2))); //将每个Bin的起始位置加入集合中
            geneticPosLists[index].add(Float.parseFloat(t.getCell(i, 3)));
        }
        List<File> fList = IOUtils.getFileListInDirEndsWith(dirS, ".txt.gz");
        fList.parallelStream().forEach(f -> {
            Dyad<String, List<String>> two = VMapDBUtils.getDBInfo(f.getAbsolutePath()); //返回Dyad类型
            String header = two.getFirstElement(); //返回表头
            List<String> recordList = two.getSecondElement(); //返回每一行的内容的集合
            String[] tem = header.split("\t");
            try {
                BufferedWriter bw = IOUtils.getTextGzipWriter(f.getAbsolutePath());
                int chrIndex = -1;
                int posIndex = -1;
                int currentPos = -1;
                String geneticPos = null;
                String ID = null;
                List<String> l  = null;
                for (int i = 0; i < recordList.size(); i++) {
                    l = PStringUtils.fastSplit(recordList.get(i)); //读每一行的内容
                    ID = l.get(0);
                    chrIndex = Integer.parseInt(l.get(1))-1; // 索引1 含有染色体号
                    currentPos = Integer.parseInt(l.get(2)); //索引2 含有位置
                    posIndex = startLists[chrIndex].binarySearch(currentPos); //在 刚刚的起始集合里搜索 index
                    if (posIndex < 0) posIndex = -posIndex-2;
                    if(posIndex == -1){ //说明在第一个数字前面
                        geneticPos = String.format("%.6f", geneticPosLists[chrIndex].get(0)/100);
                        bw.write(ID + "\t" + l.get(1)+ "\t" + geneticPos + "\t" + l.get(2)+ "\t" + l.get(3)+ "\t" + l.get(4));
                        bw.newLine();
                    }
                    else if(posIndex == startLists[chrIndex].size()-1){ //说明在最后一个数字后面
                        geneticPos = String.format("%.6f", geneticPosLists[chrIndex].get(startLists[chrIndex].size()-1)/100);
                        bw.write(ID + "\t" + l.get(1)+ "\t" + geneticPos + "\t" + l.get(2)+ "\t" + l.get(3)+ "\t" + l.get(4));
                        bw.newLine();
                    }
                    else{ //开始进行该点和集合点的比较，找最近值
                        int pos1 = currentPos - startLists[chrIndex].get(posIndex);
                        int pos2 = startLists[chrIndex].get(posIndex+1) - currentPos;
                        if(pos1 > pos2){ //说明离2近
                            geneticPos = String.format("%.6f", geneticPosLists[chrIndex].get(posIndex+1)/100);
                            bw.write(ID + "\t" + l.get(1)+ "\t" + geneticPos + "\t" + l.get(2)+ "\t" + l.get(3)+ "\t" + l.get(4));
                            bw.newLine();
                        }else{ //说明离1近
                            geneticPos = String.format("%.6f", geneticPosLists[chrIndex].get(posIndex)/100);
                            bw.write(ID + "\t" + l.get(1)+ "\t" + geneticPos + "\t" + l.get(2)+ "\t" + l.get(3)+ "\t" + l.get(4));
                            bw.newLine();

                        }

                    }
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    /**
     * 通过 Chr1A- Chr7D pos 信息转换为 chrID 格式
     * 输入文件是 genetic pos map
     *
     */
    public void convertCoordinate () {
        String infileS = "/Users/Aoyue/Documents/Data/wheat/article/iwgsc_refseqv1.0_recombination_rate_analysis/iwgsc_refseqv1.0_mapping_data.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/007_recombination/001_recombination/iwgsc_refseqv1.0_mapping_data_chrID.txt";
        RowTable<String> t = new RowTable<>(infileS);
        String header = "psId\tchromosome\tphysicalPosition\tgeneticPosition";
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(header);
            bw.newLine();
            String chromosome = null;
            String psID = null;
            String geneticPos = null;
            int pos = -1;
            int chrid = -1;
            int posid = -1;
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < t.getRowNumber(); i++) {
                sb.setLength(0);
                psID = t.getCell(i,0);
                chromosome = t.getCell(i,1).replaceFirst("chr", "");
                pos = Integer.parseInt(t.getCell(i,2));
                geneticPos = t.getCell(i,3);
                chrid = RefV1Utils.getChrID(chromosome, pos);
                posid = RefV1Utils.getPosOnChrID(chromosome, pos);
                sb.append(psID).append("\t").append(chrid).append("\t").append(posid).append("\t").append(geneticPos);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }


    /**
     * 解析老师的结果
     */
    public void extractInfoFromVMap2 () {
        int subLength = 150;
        String outDirS = "/Users/feilu/Documents/analysisH/vmap2/002_genicSNP/001_genicSNPByChr/";
        String vmapDirS = "/Volumes/Fei_HDD_Mac/VMap2.1/";
        File[] fs  = new File(vmapDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        List<File> vmapList = Arrays.asList(fs);
        Collections.sort(vmapList);
        String geneHCFileS = "/Users/feilu/Documents/analysisH/vmap2/001_geneHC/geneHC.txt";
        Table t = TablesawUtils.readTsv(geneHCFileS); //读进表格里
        System.out.println(t.structure());
        t.sortAscendingOn("Chr", "TranStart"); //升序
        IntColumn chrColumn = t.intColumn("chr"); //返回一个类 InColumn
        int chrNum = chrColumn.countUnique(); //意思是一共有42条染色体
        TIntList[] startLists = new TIntList[chrNum]; //所有的起始位点建立一个集合
        TIntList[] endLists = new TIntList[chrNum]; //所有的终止位点建立一个集合
        List<String>[] tranLists = new ArrayList[chrNum]; //每条染色体都有一个list
        for (int i = 0; i < chrNum; i++) { //对list数组进行初始化
            startLists[i] = new TIntArrayList();
            endLists[i] = new TIntArrayList();
            tranLists[i] = new ArrayList();
        }
        for (int i = 0; i < t.rowCount(); i++) {
            startLists[Integer.parseInt(t.getString(i, 2))-1].add(Integer.parseInt(t.getString(i, 3)));
            endLists[Integer.parseInt(t.getString(i, 2))-1].add(Integer.parseInt(t.getString(i, 4)));
            tranLists[Integer.parseInt(t.getString(i, 2))-1].add(t.getString(i, 1));
        }
        vmapList.parallelStream().forEach(f -> {
            int chrIndex = Integer.parseInt(f.getName().substring(3, 6))-1;
            String outfileS = new File (outDirS, f.getName().replaceFirst(".vcf.gz", "_genicSNP.txt.gz")).getAbsolutePath();
            int[] dc = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
            Arrays.sort(dc);
            StringBuilder sb = new StringBuilder();
            if (Arrays.binarySearch(dc, chrIndex+1) < 0) {
                sb.append("ID\tChr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tAAF_ABD\tAAF_AB\tTranscript");
            }
            else {
                sb.append("ID\tChr\tPos\tRef\tAlt\tMajor\tMinor\tMaf\tAAF_ABD\tAAF_D\tTranscript");
            }
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                bw.write(sb.toString());
                bw.newLine();
                String temp = null;
                while ((temp = br.readLine()).startsWith("#")) {}

                List<String> l = null;
                List<String> ll = null;
                List<String> lll = null;
                String info = null;
                int currentPos = -1;
                int posIndex = -1;
                while ((temp = br.readLine()) != null) {
                    sb.setLength(0);
                    int currentSub = subLength;
                    if (temp.length() < subLength) {
                        currentSub = temp.length();
                    }
                    l = PStringUtils.fastSplit(temp.substring(0, currentSub));
                    currentPos = Integer.parseInt(l.get(1));
                    posIndex = startLists[chrIndex].binarySearch(currentPos);
                    if (posIndex < 0) {
                        posIndex = -posIndex-2; //确保该位点在起始位点的右边
                    }
                    if (posIndex < 0) continue; //如果不在起始位点的右边，那么就不在范围内，跳过该位点
                    if (currentPos >= endLists[chrIndex].get(posIndex)) continue; //确保在末端位点的前面，若不在，也舍去
                    sb.append(l.get(2)).append("\t").append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(3));
                    sb.append("\t").append(l.get(4)).append("\t");
                    ll = PStringUtils.fastSplit(l.get(7), ";");
                    lll = PStringUtils.fastSplit(ll.get(2).replaceFirst("AD=", ""),",");
                    if (Integer.parseInt(lll.get(0)) > Integer.parseInt(lll.get(1))) {
                        sb.append(l.get(3)).append("\t").append(l.get(4)).append("\t");
                    }
                    else {
                        sb.append(l.get(4)).append("\t").append(l.get(3)).append("\t");
                    }
                    sb.append(ll.get(6).split("=")[1]).append("\t").append(ll.get(7).split("=")[1]).append("\t").append(ll.get(8).split("=")[1]);
                    sb.append("\t").append(tranLists[chrIndex].get(posIndex));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getAbsolutePath() + " is completed.");
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });

    }

}
