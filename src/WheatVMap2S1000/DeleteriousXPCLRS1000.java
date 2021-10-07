package WheatVMap2S1000;

import AoUtils.AoFile;
import gnu.trove.list.array.TCharArrayList;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.table.RowTable;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class DeleteriousXPCLRS1000 {

    public DeleteriousXPCLRS1000(){
        this.callDeleteriousXPCLR();

    }

    /**
     * 先计算不同变异类型下的个体load结果,再和同义突变进行标准化
     * @param1 choice1 选择要计算的不同变异类型的计数 choice2_ifSelected  choice3_refobj
     * @param2 parentDirS 结果文件夹
     * @param3 exonVCFDirS geneVCF vcf文件夹路径
     * @param4 SNPAnnoFileS SNP Annotation DB位点数据库，注意是42条染色体都合并在一起的文件
     * @param5 index 染色体号的索引
     * @param6 pos pos的索引号
     * @param7 variantType 数据库中的变异类型（Synonymous Nonsynonymous Noncoding Stop-loss Stop-gain Start-lost）
     * @param8 sift 即derivedSIFT值的那一列的index
     * @param9 gerp GERP信息所在的列的index
     * @param10 ancestralAllele index
     * @param11 ref index
     * @param12 alt index
     * @param13 majorAllele index
     * @param14 minorAllele index
     * @param14 snpeef vep impact index
     * @param15 定义不同变异类型的条件？？？ 也需要修改
     * @param16 vmap2TaxaList taxaInfoDB数据库
     * @param17 vcfID 中编号所在taxaInfoDB的列index
     * @param18 delVmapFileS VCF文件的名字补充
     * @param19
     *
     */
    public void callDeleteriousXPCLR(){
        //**************************** 不必修改 **************************//
        String variantType = null;
        String variantType1 = "001_synonymous";
        String variantType2 = "002_nonsynonymous";
        String variantType3 = "003_nonsynGERPandDerivedSIFT";
        String variantType4 = "004_nonsynDerivedSIFT";
        String variantType5 = "005_GERP";
        String variantType6 = "006_StopGain"; //只是SIFT中为
        String variantType7 = "007_VEP";
        String variantType8 = "008_snpEff";
        String variantType9 = "009_VEP_stopGained";
        String variantType10 = "010_GERP16way";
        String variantType11 = "011_GERP16way_1.2_max";
        String variantType12 = "012_GERP16way_1.2_2.5";
        String variantType13 = "013_GERP16way_2.5_max";
        String variantType14 = "014_GERP16way_2.14_max";
        String variantType15 = "015_LIST_S2";
        String variantType16 = "016_PhyloP";
        String variantType17 = "017_PhyloP_RefMask";
        String variantType18 = "018_GERP16wayandSIFT";
        String variantType19 = "019_Alt_PolyPhen2";
        String variantType20 = "020_Derived_PolyPhen2";
        String variantType21 = "021_GERP16wayandPolyPhen2";








        String ifselected = null;
        String ifselected1 = "1";
        String ifselected2 = "0";
        String ifselected3 = "Whole_genome";


        String group = null;
        String group1 = "wede";
        String group2 = "dedurum";
        String group3 = "lrcul";
        String group4 = "allIndivi";
        //**************************** up 不必修改 **************************//

        ////*********** selected sweep or not selected sweep *********////
//        String[] choice1_variantType = {variantType1,variantType2,variantType3,variantType4,variantType5, variantType6,variantType7,variantType8,variantType9};
//        String[] choice2_ifSelected = {ifselected1,ifselected2};
//        String[] choice3_refobj = {group1,group2,group3};

        ////*********** whole genome *********////
//        String[] choice1_variantType = {variantType1,variantType2,variantType3,variantType4,variantType5, variantType6,variantType7,variantType8,variantType9,variantType10};
//        String[] choice1_variantType = {variantType1,variantType2,variantType7,variantType11,variantType12,variantType13};
//        String[] choice1_variantType = {variantType3,variantType4,variantType5, variantType6,variantType8,variantType9,variantType10};
        String[] choice1_variantType = {variantType19,variantType20,variantType21};
        String[] choice2_ifSelected = {ifselected3};
        String[] choice3_refobj = {group4};

        ////*********** selected sweep or not selected sweep *********////
//        String[] choice1_variantType = {variantType1,variantType2,variantType3,variantType4,variantType5, variantType6,variantType7,variantType8,variantType9,variantType10};
//        String[] choice2_ifSelected = {ifselected1};
//        String[] choice3_refobj = {group1,group2,group3};



//        String parentDirS = ""; //model
//        String parentDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/009_xpclr/003_summary_XPCLR/003_deleteriousXPCLR";  //结果文件需要存放的地方
//        String parentDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/005_delCount/001_allIndivi";
//        String parentDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/005_delCount/002_allIndivi_syntenicGene";
        //从2021-08-25 开始区分 top5 top1
//        String parentDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/009_xpclr/004_summary_XPCLR_top005/003_deleteriousXPCLR";
//        String parentDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/009_xpclr/005_summary_XPCLR_top001/003_deleteriousXPCLR";

//        String parentDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/005_delCount/004_allIndivi";
        String parentDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/005_delCount/005_allIndivi";


        new File(parentDirS).mkdirs();

        for (int i = 0; i < choice3_refobj.length; i++) { //第一层循环：选择区的wede dedurum lrcl 或者是 全基因组
            group = choice3_refobj[i];
            // infileS 文件是Gene SNP 数据库中受选择区域内的 SNP 位点，即为 Annotation 数据库的子集, 如果计算的是全基因组区域的话，则该文件在程序中不会被处理到，可以忽略。
//            String infileS = new File("/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/009_xpclr/003_summary_XPCLR/002_topK","top0.05_" + group + "_ChrPos_fromExonAnnotation.txt.gz").getAbsolutePath();
//            String infileS = new File("/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/009_xpclr/004_summary_XPCLR_top005/002_topK","top0.05_" + group + "_ChrPos_fromExonAnnotation.txt.gz").getAbsolutePath();
            String infileS = new File("/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/009_xpclr/005_summary_XPCLR_top001/002_topK","top0.01_" + group + "_ChrPos_fromExonAnnotation.txt.gz").getAbsolutePath();

            for (int j = 0; j < choice2_ifSelected.length; j++) { //第二层循环：是否受选择
                ifselected = choice2_ifSelected[j];
                for (int k = 0; k < choice1_variantType.length; k++) { //第二层循环：计数的类型
                    variantType = choice1_variantType[k];
                    String DelCountFileS = new File(parentDirS,variantType + "_ifselected" + ifselected + "_" + group + "_DelCount_bychr.txt").getAbsolutePath();
                    this.countDeleteriousVMapII_byChr(infileS,ifselected,variantType,DelCountFileS,group);
                }
            }
        }
    }

    /**
     * 计算指定选择区域的mutation burden 一共有 5 步
     * @param infileS
     * @param ifselected
     * @param type
     * @param DelCountFileS
     */
    public void countDeleteriousVMapII_byChr(String infileS,String ifselected,String type,String DelCountFileS, String group) {
        //model

        //######## 需要修改 ########//
//        String exonVCFDirS = ""; //外显子变异数据
//        String SNPAnnoFileS = ""; //注释信息库合并后的总文件
//        AoFile.readheader(SNPAnnoFileS);

        int cntNONSY = 0; //非同义突变的个数

        //########### 大麦和黑麦简约法 ******* VMap2.0-2021 *********** 加上Derived SIFT的数据库 ################w
        String exonVCFDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/007_geneVCF"; //外显子变异数据 20210808 完成
//        String SNPAnnoFileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/006_geneSNPAnnotation_merge/001_geneSNPAnno.txt.gz"; //注释信息库合并后的总文件
//        String SNPAnnoFileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/006_geneSNPAnnotation_merge/002_geneSNPAnno_syntenic.txt.gz";
//        String SNPAnnoFileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/006_geneSNPAnnotation_merge/007_geneSNPAnno.txt.gz";
        String SNPAnnoFileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/006_geneSNPAnnotation_merge/011_geneSNPAnno.txt.gz";

        AoFile.readheader(SNPAnnoFileS);
/**
 * ################################### step0: 建立受选择区域的集合，并在下文进行 posList 和 ancestral charList 构建时进行适当的过滤。
 */

        TIntArrayList[] selectedPosList = new TIntArrayList[42]; //不同染色体下，不同受选择区域位点的集合

        if (!ifselected.equals("Whole_genome")){
            // I only need to know the chr pos information about the selected region
            for (int i = 0; i < selectedPosList.length; i++) {
                selectedPosList[i] = new TIntArrayList();
            }
            RowTable<String> selectedT = new RowTable<>(infileS);
            for (int i = 0; i < selectedT.getRowNumber(); i++) {
                int chr = selectedT.getCellAsInteger(i,0);
                int pos = selectedT.getCellAsInteger(i,1);
                selectedPosList[chr-1].add(pos);
            }
            for (int i = 0; i < selectedPosList.length; i++) {
                selectedPosList[i].sort();
            }
            System.out.println("Finished step0: completing the initialization of selected list.");
        }

        /**
         *  ################################### step1: 初始化染色体集合
         */
        int minDepth = 2;//inclusive
        int chrNum = 42;
        ArrayList<Integer> chrList = new ArrayList();
        for (int i = 0; i < chrNum; i++) {
            chrList.add(i + 1);
        }
        System.out.println("Finished step1: completing the initialization of chromosome.");

        /**
         *  ################################### step2: posList  charList 补充完整
         */

        int[][] delePos = new int[chrNum][];
        char[][] deleChar = new char[chrNum][];

        TIntArrayList[] posList = new TIntArrayList[chrNum];
        TCharArrayList[] charList = new TCharArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) { //集合类数组，要初始化每一个list
            posList[i] = new TIntArrayList();
            charList[i] = new TCharArrayList();
        }

        String derivedAllele = null;

        try {
            BufferedReader br = AoFile.readFile(SNPAnnoFileS);
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                int index = Integer.parseInt(l.get(1)) - 1; //染色体号的索引 ################ 需要修改 需要修改 需要修改 ################
                int pos = Integer.parseInt(l.get(2)); //################ 需要修改 需要修改 需要修改 ################

                /**
                 * 进行受选择区域的判断,如果受选择，那么判断该位点是否受选择，不受选择的进行过滤。
                 * 如果不受选择，那么判断该位点是否受选择，受选择的进行过滤。
                 */

                if (!ifselected.equals("Whole_genome")){
                    int indexselect = selectedPosList[index].binarySearch(pos);
                    if (ifselected.equals("1")){ //只保留受选择的位点，把不受选择的进行过滤
                        if (indexselect <0) continue;
                    }
                    if (ifselected.equals("0")){ //只保留不受选择的位点，把受选择的进行过滤
                        if (indexselect>= 0) continue;
                    }
                }


                String variantType = l.get(13); //################ 需要修改 需要修改 需要修改 ################
                String sift = l.get(16); //################ 需要修改 需要修改 需要修改 ################
                String gerp = l.get(11); //################ 需要修改 需要修改 需要修改 ################
                String gerp16way = l.get(22);

                //********************* 过滤没有 ancestral allele 信息的位点
                //################### 需要修改 //###################//###################//###################//###################
                ////不同的数据库，这一列的信息不一样，千万要注意!!!!!!!!!!!!!!!!! 祖先基因的数据库
                String ancestralAllele = l.get(9);
                //################### 需要修改 //###################//###################//###################//###################
                String ref = l.get(3);
                String alt = l.get(4);
                String majorAllele = l.get(5);
                String minorAllele = l.get(6);
                String Effect_VEP = l.get(17);
                String Impact_VEP = l.get(18);
                String Impact_snpEff = l.get(20);
                String LIST_S2 = l.get(25);
                String phyloP = l.get(26);
                String phyloP_RefMask = l.get(28);
                String altPolyPhen = l.get(34);
                String derivedPolyPhen = l.get(40);
                if(!(ancestralAllele.equals(ref)||ancestralAllele.equals(alt)))continue;


                /**
                 ******** 定义有害突变，不是有害突变，就忽略不计 ################ 需要修改 需要修改 需要修改 ################
                 */

                if (type.equals("001_synonymous")){
                    if (!variantType.equals("SYNONYMOUS"))continue;
                }
                if (type.equals("002_nonsynonymous")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue;
                }
                if (type.equals("003_nonsynGERPandDerivedSIFT")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                    if(gerp.startsWith("N")) continue; //说明必须满足GERP有值
                    if (sift.startsWith("N"))continue; //说明必须满足有sift值
                    double gerpd = Double.parseDouble(gerp);
                    double siftd = Double.parseDouble(sift);
                    if (gerpd <= 1) continue; //说明必须满足gerp大于1
                    if (siftd >= 0.05 ) continue; //说明必须满足sift小于0.05
                }

                if (type.equals("004_nonsynDerivedSIFT")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                    if (sift.startsWith("N"))continue; //说明必须满足有sift值
                    double siftd = Double.parseDouble(sift);
                    if (siftd >= 0.05) continue; //说明必须满足sift小于0.05
                }

                if (type.equals("005_GERP")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                    if(gerp.startsWith("N")) continue; //说明必须满足GERP有值
                    double gerpd = Double.parseDouble(gerp);
                    if (gerpd <= 1) continue; //说明必须满足gerp大于1
                }

                if(type.equals("006_StopGain")){
                    if (!variantType.equals("STOP-GAIN") && !variantType.equals("STOP-LOSS") && !variantType.equals("START-LOST"))continue; //等于三者之一
                }

                if(type.equals("007_VEP")){
                    if (!Impact_VEP.equals("HIGH"))continue; //
                }

                if(type.equals("009_VEP_stopGained")){
                    if (!Effect_VEP.contains("start_lost") && !Effect_VEP.contains("stop_gained") && !Effect_VEP.contains("stop_lost")) continue;
                }

                if(type.equals("008_snpEff")){
                    if (!Impact_snpEff.equals("HIGH"))continue; //
                }

                if(type.equals("010_GERP16way")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                    if(gerp16way.startsWith("N")) continue; //说明必须满足GERP有值
                    double gerp16wayd = Double.parseDouble(gerp16way);
                    if (gerp16wayd <= 1) continue; //说明必须满足gerp大于等于1
                }

                if(type.equals("011_GERP16way_1.2_max")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                    if(gerp16way.startsWith("N")) continue; //说明必须满足GERP有值
                    double gerp16wayd = Double.parseDouble(gerp16way);
                    if (gerp16wayd <= 1.2) continue; //说明必须满足gerp大于等于1
                }

                if(type.equals("012_GERP16way_1.2_2.5")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                    if(gerp16way.startsWith("N")) continue; //说明必须满足GERP有值
                    double gerp16wayd = Double.parseDouble(gerp16way);
                    if (gerp16wayd <= 1.2 || gerp16wayd >= 2.5) continue; //说明必须满足gerp大于等于1
                }

                if(type.equals("013_GERP16way_2.5_max")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                    if(gerp16way.startsWith("N")) continue; //说明必须满足GERP有值
                    double gerp16wayd = Double.parseDouble(gerp16way);
                    if (gerp16wayd <= 2.5) continue; //说明必须满足gerp大于等于1
                }

                if(type.equals("014_GERP16way_2.14_max")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                    if(gerp16way.startsWith("N")) continue; //说明必须满足GERP有值
                    double gerp16wayd = Double.parseDouble(gerp16way);
                    if (gerp16wayd < 2.14) continue; //说明必须满足gerp大于等于1
                }

                if(type.equals("015_LIST_S2")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                    if(LIST_S2.startsWith("N")) continue; //说明必须满足 LIST 有值
                    double LIST_S2d = Double.parseDouble(LIST_S2);
                    if (LIST_S2d < 0.85) continue; //说明必须满足 list 大于等于 0.85
                }

                if(type.equals("016_PhyloP")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                    if(phyloP.startsWith("N")) continue; //说明必须满足有值
                    double phyloPd = Double.parseDouble(phyloP);
                    if (phyloPd < 1.5) continue; //说明必须满足 list 大于等于 0.85
                }

                if(type.equals("017_PhyloP_RefMask")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                    if(phyloP_RefMask.startsWith("N")) continue; //说明必须满足有值
                    double phyloP_RefMaskd = Double.parseDouble(phyloP_RefMask);
                    if (phyloP_RefMaskd < 1.5) continue; //说明必须满足 list 大于等于 0.85
                }

                if (type.equals("018_GERP16wayandSIFT")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                    if(gerp16way.startsWith("N")) continue; //说明必须满足GERP有值
                    if (sift.startsWith("N"))continue; //说明必须满足有sift值
                    double gerp16wayd = Double.parseDouble(gerp16way);
                    double siftd = Double.parseDouble(sift);
                    if (gerp16wayd < 2.14) continue; //说明必须满足gerp大于等于2.14
                    if (siftd >= 0.05 ) continue; //说明必须满足sift小于0.05
                }

                if (type.equals("019_Alt_PolyPhen2")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                    if (altPolyPhen.startsWith("N"))continue; //说明必须满足有sift值
                    double altPolyPhend = Double.parseDouble(altPolyPhen);
                    if (altPolyPhend < 0.453) continue; //说明必须满足gerp大于等于2.14
                }

                if (type.equals("020_Derived_PolyPhen2")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                    if (derivedPolyPhen.startsWith("N"))continue; //说明必须满足有sift值
                    double derivedPolyPhend = Double.parseDouble(derivedPolyPhen);
                    if (derivedPolyPhend < 0.453) continue; //说明必须满足gerp大于等于2.14
                }

                if (type.equals("021_GERP16wayandPolyPhen2")){
                    if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                    if(gerp16way.startsWith("N")) continue; //说明必须满足GERP有值
                    if (derivedPolyPhen.startsWith("N"))continue; //说明必须满足有sift值
                    double gerp16wayd = Double.parseDouble(gerp16way);
                    double derivedPolyPhend = Double.parseDouble(derivedPolyPhen);
                    if (gerp16wayd < 2.14) continue; //说明必须满足gerp大于等于2.14
                    if (derivedPolyPhend < 0.453 ) continue; //说明必须满足sift小于0.05
                }

                //String variantType18 = "018_GERP16wayandSIFT";

                //        String variantType17 = "017_PhyloP_RefMask";

                //        String variantType15 = "015_LIST_S2";
                //        String variantType16 = "016_PhyloP";

                ////////////////////// depracate
                if (type.equals("006_nonsynGERPandDerivedSIFT_correction")){
                    if (ancestralAllele.equals(ref)){ //
                        if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                        if(gerp.startsWith("N")) continue; //说明必须满足GERP有值
                        if (sift.startsWith("N"))continue; //说明必须满足有sift值
                        double gerpd = Double.parseDouble(gerp);
                        double siftd = Double.parseDouble(sift);
                        if (gerpd < 1.61) continue; //说明必须满足gerp大于1
                        if (siftd > 0.04) continue; //说明必须满足sift小于等于0.05
                    }
                    if (ancestralAllele.equals(alt)){
                        if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                        if(gerp.startsWith("N")) continue; //说明必须满足GERP有值
                        if (sift.startsWith("N"))continue; //说明必须满足有sift值
                        double gerpd = Double.parseDouble(gerp);
                        double siftd = Double.parseDouble(sift);
                        if (gerpd < 0.0139) continue; //说明必须满足gerp大于1
                        if (siftd > 0.3) continue; //说明必须满足sift小于等于0.05
                    }
                }

                if (type.equals("007_nonsynDerivedSIFT_correction")){
                    if (ancestralAllele.equals(ref)){ //
                        if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                        if (sift.startsWith("N"))continue; //说明必须满足有sift值
                        double siftd = Double.parseDouble(sift);
                        if (siftd > 0.04) continue; //说明必须满足sift小于等于0.05
                    }
                    if (ancestralAllele.equals(alt)){
                        if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                        if (sift.startsWith("N"))continue; //说明必须满足有sift值
                        double siftd = Double.parseDouble(sift);
                        if (siftd > 0.3) continue; //说明必须满足sift小于等于0.05
                    }
                }

                if (type.equals("008_GERP_correction")){
                    if (ancestralAllele.equals(ref)){ //
                        if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                        if(gerp.startsWith("N")) continue; //说明必须满足GERP有值
                        double gerpd = Double.parseDouble(gerp);
                        if (gerpd < 1.61) continue; //说明必须满足gerp大于1
                    }
                    if (ancestralAllele.equals(alt)){
                        if (!variantType.equals("NONSYNONYMOUS"))continue; //说明必须满足是非同义突变
                        if(gerp.startsWith("N")) continue; //说明必须满足GERP有值
                        double gerpd = Double.parseDouble(gerp);
                        if (gerpd < 0.0139) continue; //说明必须满足gerp大于1
                    }
                } // end


                if (ancestralAllele.equals(majorAllele)) {
                    derivedAllele = minorAllele;
                    posList[index].add(pos); //将包含有derived allele的位点添加到Poslist
                    charList[index].add(derivedAllele.charAt(0)); //返回的是char类型的字符
                }
                if (ancestralAllele.equals(minorAllele)) {
                    derivedAllele = majorAllele;
                    posList[index].add(pos);
                    charList[index].add(derivedAllele.charAt(0));
                }
                else if (!(ancestralAllele.equals(majorAllele) || ancestralAllele.equals(minorAllele))){
                }
                cnt++;
                if (cnt%100000==0) System.out.println("cnt is " + cnt +" going on at step2");

            }
            br.close();
            System.out.println(cntNONSY + "SNP num in " + type);

            for (int i = 0; i < chrNum; i++) { //将每一个list转化为数组
                delePos[i] = posList[i].toArray();
                deleChar[i] = charList[i].toArray();
                Arrays.sort(delePos[i]);
            }
            System.out.println("Finished step2: completing the posList  charList.");

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        /**
         *  ################################### step3: taxa 集合 642个taxa
         */
        String vmap2TaxaList = "/Users/Aoyue/project/wheatVMap2_1000/001_germplasm/021_WheatVMap2_GermplasmInfo.txt";

        String[] taxa = AoFile.getStringArraybyList(vmap2TaxaList,0);

        /**
         * count 计数有三种类型：1.total, 2.纯合子， 3.杂合子
         */
        double[][] totalDelCount = new double[chrNum][taxa.length];
        int[][] homoDelCount = new int[chrNum][taxa.length];
        double[][] heterDelCount = new double[chrNum][taxa.length];
        int[][] siteWithMinDepthCount = new int[chrNum][taxa.length]; //每个taxa在每条染色体中的有害突变位点

        chrList.parallelStream().forEach(chr -> {
            String delVmapFileS = "chr" + PStringUtils.getNDigitNumber(3, chr) + "_gene_vmap2.1.vcf.gz";
            //开始读写VCF文件
            delVmapFileS = new File(exonVCFDirS, delVmapFileS).getAbsolutePath();
            BufferedReader br = AoFile.readFile(delVmapFileS);
            int chrIndex = chr - 1;
            try {
                String temp = null;
                int cnt = 0;
                HashMap<String, Integer> hmtaxainVCFindex = new HashMap<>();
                HashMap<String, Integer> hmtaxainTaxaindex = new HashMap<>();
                //这里涉及2个数组概念， taxa and taxainVCFfile,即总的taxa和在VCF文件中的taxa
                String[] taxainVCFfile = null;
                List<String> taxainVCFlist = new ArrayList();
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("##")) continue;
                    if (temp.startsWith("#CHR")) {//说明进入taxa信息列
                        List<String> l = PStringUtils.fastSplit(temp, "\t");
                        taxainVCFfile = new String[l.size() - 9];
                        for (int i = 9; i < l.size(); i++) {
                            taxainVCFfile[i - 9] = l.get(i);
                            hmtaxainVCFindex.put(taxainVCFfile[i - 9], i); //第0个taxa的genotype在第9列，第一个taxa的genotype在第10列，依次类推；
                        }
                        // ************** 在总的taxa中，搜索VCF文件中的taxa的index
                        for (int i = 0; i < taxainVCFfile.length; i++) { //找到taxa
                            int index = Arrays.binarySearch(taxa, taxainVCFfile[i]);
                            if (index > -1) {
                                hmtaxainTaxaindex.put(taxainVCFfile[i], index);
                            }
                        }
                        taxainVCFlist = Arrays.asList(taxainVCFfile);
                    }

                    if (!temp.startsWith("#")) {
                        cnt++;
                        if (cnt % 10000 == 0) {
                            System.out.println(String.valueOf(cnt/1000) + " kb lines on chr " + String.valueOf(chr));
                        }
                        List<String> l = PStringUtils.fastSplit(temp.substring(0, 100), "\t");
                        int pos = Integer.valueOf(l.get(1)); //根据pos找到index
                        int index = Arrays.binarySearch(delePos[chrIndex], pos);
                        if (index < 0) {
                            continue; //前面建立的 delePos[][] 和deleChar[][] 都是为现在在vcf文件中找位置贡献的，不是有害突变的位点，都过滤。*************************************************************
                        } //重要：假如该位点有基因型，也是非同义突变位点，但是没有祖先位点状态，于是就过滤。
                        l = PStringUtils.fastSplit(temp, "\t");
                        int[] idx = new int[2];
                        if (l.get(3).charAt(0) == deleChar[chrIndex][index]) { //如果ref allele = deleChar allele
                            idx[0] = 0; //等于0则代表是 derived allele
                            idx[1] = 1;
                        } else {
                            idx[0] = 1;
                            idx[1] = 0;
                        }

                        //合计642个taxa，在A Bgenome中只有606（814+212=1026）个，在Dgenome中只有814+36=850个，我们要找到每个VCF文件中的genotype所对应的taxa的index
                        //code:本段代码是每行SNP位点，每个taxa的有害位点的统计
                        for (int i = 0; i < taxainVCFlist.size(); i++) {
                            int genotypeIndex = hmtaxainVCFindex.get(taxainVCFlist.get(i)); //获取该taxa的基因型所在的列的索引
                            int taxaIndex = hmtaxainTaxaindex.get(taxainVCFlist.get(i)); //获取该taxa所在总的642个数组中的索引，为后续写文件进行统计
                            String genoS = l.get(genotypeIndex); //指的是 GT:AD:GL 信息
                            if (genoS.startsWith(".")) {
                                continue; //如果以.开头，说明没有基因型信息，此位点没有测到。
                            }
                            List<String> ll = PStringUtils.fastSplit(genoS, ":"); //分开为GT   AD   GL三类

                            List<String> lll = PStringUtils.fastSplit(ll.get(1), ","); //lll指将AD提取出来，并以"，"号分割。如 0/0:1,0:0,3,28中 ，1，0分别代表ref和alt的测序深度
                            int depth = Integer.valueOf(lll.get(0)) + Integer.valueOf(lll.get(1)); //总得测序深度等于 ref + alt
                            if (depth < minDepth) continue; //最小测序深度是2，如果小于2，则弃用
                            lll = PStringUtils.fastSplit(ll.get(0), "/"); //这里lll指的是基因型GT，lll被重新赋值，之前代表的是AD
                            int v1 = Integer.valueOf(lll.get(0)); //v1等于 ref
                            int v2 = Integer.valueOf(lll.get(1)); // v2 等于 alt
                            int sum = 0;
                            //如果ref是derived allele，那么idx[0]=0;idx[1]=1. 当是0/0时，sum=2; 0/1时，sum=1; 1/1时， sum=0.
                            //如果alt是derived allele，那么idx[0]=1;idx[1]=0. 当是0/0时，sum=0; 0/1时，sum=1; 1/1时， sum=2.
                            //      巧妙利用 0/0 和 derived allele 的关系，如果 ref 是 derived, 那么 idx[0] = 0， 效应可以加1; 如果 ref 是 ancestral, 那么 idx[0] = 1, 效应就不能加1了;
                            //
                            if (v1 == idx[0]) { //
                                sum++;
                            }
                            if (v2 == idx[0]) {
                                sum++;
                            }
                            if (sum == 0) {
                            } else if (sum == 1) {
                                totalDelCount[chrIndex][taxaIndex] += 0.5;
                                heterDelCount[chrIndex][taxaIndex] += 0.5;
                            } else {
                                totalDelCount[chrIndex][taxaIndex] += 1;
                                homoDelCount[chrIndex][taxaIndex] += 1;
                            }
                            siteWithMinDepthCount[chrIndex][taxaIndex]++; //taxa有多少个有害突变位点
                        }
                    }
                }
                br.close();
                System.out.println("Finished step3: calculating each taxa for each chromosome.");

            } catch (Exception e) {
                e.printStackTrace();
            }

            /**
             * 进行每个样品每条染色体的计算
             */
        });

        /**
         *  ################################### step4: 输出每个样品每条染色体的结果
         */

        try {
            BufferedWriter bw = AoFile.writeFile(DelCountFileS);
            bw.write("Taxa\tChr\tTotalDerivedSNPCount\tHomoDerivedSNPCount\tHeterDerivedSNPCount\tSiteCountWithMinDepth\tVariantsGroup\tIfSelected\tGroup_refobj"); //每个taxa有多少个加性效应的derivedAllele, 每个taxa在del库中含有基因型的个数
            bw.newLine();
            for (int i = 0; i < homoDelCount.length; i++) { //第一层是染色体号
                String chr = String.valueOf(i + 1);
                for (int j = 0; j < homoDelCount[0].length; j++) { //第二层是该号染色体的有害变异个数
                    bw.write(taxa[j] + "\t" + chr + "\t" + String.valueOf(totalDelCount[i][j])+ "\t" + String.valueOf(homoDelCount[i][j]) + "\t" + String.valueOf(heterDelCount[i][j])+ "\t" + String.valueOf(siteWithMinDepthCount[i][j])
                            + "\t" + type + "\t" + ifselected+ "\t" + group);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();

            System.out.println("Finished step4: writing taxa table.");

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}