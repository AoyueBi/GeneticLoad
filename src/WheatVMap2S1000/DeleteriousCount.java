package WheatVMap2S1000;

import AoUtils.AoFile;
import gnu.trove.list.array.TCharArrayList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.table.RowTable;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

public class DeleteriousCount {

    /**
     * 2021年7月中旬，对VMap2.0样本增加到1000+,此时重新进行个体deleterious load的计算
     */
    public DeleteriousCount(){
//        AoFile.readheader("/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/005_SNP_Annotation_filterN/chr036_gene_geneSNP.txt.gz");
//        AoFile.readheader("/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/006_geneSNPAnnotation_merge/002_geneSiteAnno.txt.gz");
        this.callDeleterious();
//        AoFile.readheader("/Users/Aoyue/project/wheatVMap2_1000/001_germplasm/009_WheatVMap2_GermplasmInfo_20210708.txt");

    }

    /**
     * 先计算不同变异类型下的个体load结果,再和同义突变进行标准化
     * @param1 choice1 选择要计算的不同变异类型的计数
     * @param2 parentDirS 结果文件夹
     * @param3 exonVCFDirS geneVCF 文件夹路径
     * @param4 SNPAnnoFileS SNP位点数据库，注意是42条染色体都合并在一起的文件
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
    public void callDeleterious(){
//**************************** 不必修改 **************************//
        String variantType = null;
        String variantType1 = "001_synonymous";
        String variantType2 = "002_nonsynonymous";
        String variantType3 = "003_nonsynGERPandDerivedSIFT";
        String variantType4 = "004_nonsynDerivedSIFT";
        String variantType5 = "005_GERP";
        String variantType6 = "006_StopGain";
        String variantType7 = "007_VEP";
        String variantType8 = "008_snpEff";


//        String variantType7 = "006_nonsynGERPandDerivedSIFT_correction";
//        String variantType8 = "007_nonsynDerivedSIFT_correction";
//        String variantType9 = "008_GERP_correction";

        String ratioType = null;
        String ratioType1 = "bySub";
        String ratioType2 = "bysub_mergeByTaxa";
        String group = null;


        String[] choice1 = {variantType1,variantType2,variantType3,variantType4,variantType5,variantType6,variantType7,variantType8};
//        String[] choice1 = {variantType6};

        String[] choice3 = {ratioType1,ratioType2};
        //**************************** up 不必修改 **************************//


          //********** 这一步是进行 不同类型的load 的计算 ************//

//        String parentDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/005_delCount/001_allIndivi"; //结果文件需要存放的地方
        String parentDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/005_delCount/002_allIndivi_syntenicGene";
        new File(parentDirS).mkdirs();
        for (int k = 0; k < choice1.length; k++) { //每种变异类型，都输出了一个delCount的文件，到时候把文件
            variantType = choice1[k];
            String DelCountFileS = new File(parentDirS,variantType + "_DelCount_bychr.txt").getAbsolutePath();
            this.countDeleteriousVMapII_byChr(variantType,DelCountFileS);
        }

//        //********** 这一步是进行 ratio = nonsyn / syn 的计算 ************//
//        int cnt = 0;
//        for (int k = 1; k < choice1.length; k++) {
//            variantType = choice1[k];
//            for (int l = 0; l < choice3.length; l++) {
//                ratioType = choice3[l];
//                String infileS1 = new File(parentDirS,variantType + "_bychr_" + ratioType + ".txt").getAbsolutePath();
//                String infileS2 = new File(parentDirS,"001_synonymous" + "_bychr_" + ratioType +".txt").getAbsolutePath();
//                this.DeltoSynonymousRatio(infileS1,ratioType,infileS2);
//                cnt++;
//                System.out.println("********" + cnt + " " + variantType + " " + ratioType);
//            }
//        }

    }

    public void countDeleteriousVMapII_byChr(String type, String DelCountFileS) {
        //model

        //######## 需要修改 ########//
//        String exonVCFDirS = ""; //外显子变异数据
//        String SNPAnnoFileS = ""; //注释信息库合并后的总文件
//        AoFile.readheader(SNPAnnoFileS);

        int cntNONSY = 0; //非同义突变的个数

        //########### 大麦和黑麦简约法 ******* VMap2.0-2020 *********** 加上Derived SIFT的数据库 2020-07-21 ################w
        String exonVCFDirS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/007_geneVCF"; //外显子变异数据
//        String SNPAnnoFileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/006_geneSNPAnnotation_merge/001_geneSNP_anno.txt.gz"; //注释信息库合并后的总文件
//        String SNPAnnoFileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/006_geneSNPAnnotation_merge/002_geneSiteAnno.txt.gz";
        String SNPAnnoFileS = "/Users/Aoyue/project/wheatVMap2_1000/002_dataAnalysis/004_annoDB/006_geneSNPAnnotation_merge/003_geneSiteAnno_syntenic.txt.gz";

        AoFile.readheader(SNPAnnoFileS);
        //************* 无需修改的路径 ****************** //
        String DelCountFiletempS = new File(DelCountFileS).getAbsolutePath().replaceFirst(".txt",".temp.txt"); //有害变异加性模型输出文件


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
                String variantType = l.get(13); //################ 需要修改 需要修改 需要修改 ################
                String sift = l.get(16); //################ 需要修改 需要修改 需要修改 ################
                String gerp = l.get(11); //################ 需要修改 需要修改 需要修改 ################

                //********************* 过滤没有 ancestral allele 信息的位点
                //################### 需要修改 //###################//###################//###################//###################
                ////不同的数据库，这一列的信息不一样，千万要注意!!!!!!!!!!!!!!!!! 祖先基因的数据库
                String ancestralAllele = l.get(9);
                //################### 需要修改 //###################//###################//###################//###################
                String ref = l.get(3);
                String alt = l.get(4);
                String majorAllele = l.get(5);
                String minorAllele = l.get(6);
                String Impact_VEP = l.get(18);
                String Impact_snpEff = l.get(20);
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
                    if (gerpd <= 1) continue; //说明必须满足gerp大于等于1
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
                    if (gerpd <= 1) continue; //说明必须满足gerp大于等于1
                }

                if(type.equals("006_StopGain")){
                    if (!variantType.equals("STOP-GAIN") && !variantType.equals("STOP-LOSS") && !variantType.equals("START-LOST"))continue; //等于三者之一
                }

                if(type.equals("007_VEP")){
                    if (!Impact_VEP.equals("HIGH"))continue; //
                }

                if(type.equals("008_snpEff")){
                    if (!Impact_snpEff.equals("HIGH"))continue; //
                }



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
                }


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
        String vmap2TaxaList = "/Users/Aoyue/project/wheatVMap2_1000/001_germplasm/009_WheatVMap2_GermplasmInfo_20210708.txt";

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
                            hmtaxainVCFindex.put(taxainVCFfile[i - 9], i); //第0个taxa的genotype在第9行，第一个taxa的genotype在第10行，依次类推；
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


                        //合计642个taxa，在A Bgenome中只有606（419+187）个，在Dgenome中只有455（419+36）个，我们要找到每个VCF文件中的genotype所对应的taxa的index
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
                            //
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
            bw.write("Taxa\tChr\tTotalDerivedSNPCount\tHomoDerivedSNPCount\tHeterDerivedSNPCount\tSiteCountWithMinDepth\tVariantsGroup"); //每个taxa有多少个加性效应的derivedAllele, 每个taxa在del库中含有基因型的个数
            bw.newLine();
            for (int i = 0; i < homoDelCount.length; i++) { //第一层是染色体号
                String chr = String.valueOf(i + 1);
                for (int j = 0; j < homoDelCount[0].length; j++) { //第二层是该号染色体的有害变异个数
                    bw.write(taxa[j] + "\t" + chr + "\t" + String.valueOf(totalDelCount[i][j])+ "\t" + String.valueOf(homoDelCount[i][j]) + "\t" + String.valueOf(heterDelCount[i][j])+ "\t" + String.valueOf(siteWithMinDepthCount[i][j])
                            + "\t" + type);
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

    public void DeltoSynonymousRatio(String infileS1, String ratioType, String infileS2){
//        String infileS1 = ""; //有害突变文件
//        String infileS2 = ""; //同义突变文件

        String outfileS = new File(infileS1).getAbsolutePath().replaceFirst(".txt","_delVSsynonymous.txt");
        AoFile.readheader(infileS1);

        TDoubleArrayList del = new TDoubleArrayList();
        TDoubleArrayList syn = new TDoubleArrayList();
        TDoubleArrayList delcount = new TDoubleArrayList();
        TDoubleArrayList syncount = new TDoubleArrayList();

        if(ratioType.equals("bySub")){
            del = AoFile.getTDoubleList(infileS1,7); //bySub的情况
            syn = AoFile.getTDoubleList(infileS2,7); //bySub的情况
            delcount = AoFile.getTDoubleList(infileS1,2); //bySub的情况
            syncount = AoFile.getTDoubleList(infileS2,2); //bySub的情况


        }else if (ratioType.equals("bysub_mergeByTaxa")){
            del = AoFile.getTDoubleList(infileS1,6); //byTaxa的情况
            syn = AoFile.getTDoubleList(infileS2,6); //byTaxa的情况
            delcount = AoFile.getTDoubleList(infileS1,1); //byTaxa的情况
            syncount = AoFile.getTDoubleList(infileS2,1); //byTaxa的情况
        }

        TDoubleArrayList ratioList = new TDoubleArrayList(del.size());
        for (int i = 0; i < del.size(); i++) {
            double r = (double)del.get(i)/syn.get(i); //String.format("%.3f",r)
            ratioList.add(r);
        }

        //用count除以count
        TDoubleArrayList ratioList2 = new TDoubleArrayList(delcount.size());
        for (int i = 0; i < delcount.size(); i++) {
            double r = (double)delcount.get(i)/syncount.get(i); //String.format("%.3f",r)
            ratioList2.add(r);
        }

        try {
            BufferedReader br = AoFile.readFile(infileS1);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = br.readLine();
            bw.write(temp + "\tRatio_delVSsyn\tRatio_delVSsyn_bycount");
            bw.newLine();
            int i = 0;
            while ((temp = br.readLine()) != null) {
                bw.write(temp + "\t" + String.format("%.3f",ratioList.get(i))+ "\t" + String.format("%.3f",ratioList2.get(i)));
                bw.newLine();
                i++;
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

}
