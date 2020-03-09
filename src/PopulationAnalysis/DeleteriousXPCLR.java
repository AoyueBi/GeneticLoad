package PopulationAnalysis;

import AoUtils.AoFile;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TCharArrayList;
import gnu.trove.list.array.TIntArrayList;
import pgl.format.table.RowTable;
import pgl.utils.IOUtils;
import pgl.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

public class DeleteriousXPCLR {

    public DeleteriousXPCLR(){
//        this.mergeExonSNPAnnotation();
        this.countDeleteriousVMapII_byChr();
//        this.test();

    }

    public void test(){
        AoFile.readheader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/taxaList.txt");


    }

    /**
     * 计算受选择区域的 Clutivar 和 landrace Europe 之间的 mutation burden
     */
    public void countDeleteriousVMapII_byChr() {
        String delVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/002_exonSNPVCF"; //有害变异的VCF文件路径
        String deleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/004_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz"; //有害变异信息库

        //Top K xpclr regions
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/005_out/001_CLvsLR/004_merge/001_CLvsEU_exonRegion_0.0001_200_50000_addHeader_sortbyXPCLR_top0.05.xpclr.txt";
        // 受选择区域的位点列表
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/008_deleteriousRegion/001_selectedRegion/001_ExonSNP_anno_selectedRegion.txt";

//        String addCountFileS = ""; //有害变异加性模型输出文件
//        String recCountFileS = ""; //有害变异隐形模型输出文件

        String addCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/008_deleteriousRegion/002_countDel/001_additiveDeleterious_vmap2_bychr_selectedRegion.txt";
        String recCountFileAddGroupS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/008_deleteriousRegion/002_countDel/002_recessiveDeleterious_vmap2_bychr_selectedRegion.txt";

        String addCountFileS = new File(addCountFileAddGroupS).getAbsolutePath().replaceFirst(".txt",".temp.txt"); //有害变异加性模型输出文件
        String recCountFileS = new File(recCountFileAddGroupS).getAbsolutePath().replaceFirst(".txt",".temp.txt"); //有害变异隐形模型输出文件

        AoFile.readheader(deleFileS);

        /**
         *  ################################### 初始化染色体集合
         */
        int minDepth = 2;//inclusive
        int chrNum = 42;
        ArrayList<Integer> chrList = new ArrayList();
        for (int i = 0; i < chrNum; i++) {
            chrList.add(i + 1);
        }

        /**
         * 建立 region 的集合
         */
        int bin = 50000;
        int currentPos = -1;
        int chrInd = -1;
        RowTable<String> topt = new RowTable<>(infileS);
        TIntList[] startLists = new TIntList[chrNum]; //所有的起始位点建立一个集合
        TIntList[] endLists = new TIntList[chrNum]; //所有的终止位点建立一个集合
        for (int i = 0; i < chrNum; i++) { //对list数组进行初始化
            startLists[i] = new TIntArrayList();
            endLists[i] = new TIntArrayList();
        }
        for (int i = 0; i < topt.getRowNumber(); i++) {
            chrInd = topt.getCellAsInteger(i,0) -1;
            currentPos = topt.getCellAsInteger(i,3);
            startLists[chrInd].add(currentPos);
            endLists[chrInd].add(currentPos+bin);
        }
        for (int i = 0; i < startLists.length; i++) {
            startLists[i].sort();
            endLists[i].sort();
        }
        System.out.println("Finished building the region list");


        /**
         *  ################################### posList  charList 补充完整
         */

        int[][] delePos = new int[chrNum][];
        char[][] deleChar = new char[chrNum][];

        TIntArrayList[] posList = new TIntArrayList[chrNum];
        TCharArrayList[] charList = new TCharArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) { //集合类数组，要初始化每一个list
            posList[i] = new TIntArrayList();
            charList[i] = new TCharArrayList();
        }
        new AoFile().readheader(deleFileS);
        String derivedAllele = null;

        try{
            RowTable<String> t = new RowTable(deleFileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            bw.write("CHROM\tPOS");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) {
                int index = t.getCellAsInteger(i, 1) - 1; //染色体号的索引
                int pos = t.getCellAsInteger(i,2);
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
                bw.write(index+1 + "\t" + pos);
                bw.newLine();

                /**
                 * 定义有害突变，不是有害突变，就忽略不计
                 */
                String variantType = t.getCell(i,12);
                String sift = t.getCell(i,13);
                String gerp = t.getCell(i,20);
                if (!variantType.equals("NONSYNONYMOUS") || sift.equals("NA") || gerp.equals("NA"))continue;
                double siftd = Double.parseDouble(sift);
                double gerpd = Double.parseDouble(gerp);
                if (siftd >= 0.05 || gerpd <= 1)continue;
                //################### 需要修改 //###################//###################//###################//###################
//            String ancestralAllele = t.getCell(i, 22); //不同的数据库，这一列的信息不一样，千万要注意!!!!!!!!!!!!!!!!! 祖先基因的数据库
                String ancestralAllele = t.getCell(i, 15); //不同的数据库，这一列的信息不一样，千万要注意!!!!!!!!!!!!!!!!! 祖先基因的数据库
                //################### 需要修改 //###################//###################//###################//###################

                String majorAllele = t.getCell(i, 5);
                String minorAllele = t.getCell(i, 6);
                if (ancestralAllele.equals(majorAllele)) {
                    derivedAllele = minorAllele;
                    posList[index].add(t.getCellAsInteger(i, 2)); //将包含有derived allele的位点添加到Poslist
                    charList[index].add(derivedAllele.charAt(0)); //返回的是char类型的字符
                }
                if (ancestralAllele.equals(minorAllele)) {
                    derivedAllele = majorAllele;
                    posList[index].add(t.getCellAsInteger(i, 2));
                    charList[index].add(derivedAllele.charAt(0));
                }
                else if (!(ancestralAllele.equals(majorAllele) || ancestralAllele.equals(minorAllele))){
                }
            }
            bw.flush();
            bw.close();

            for (int i = 0; i < chrNum; i++) { //将每一个list转化为数组
                delePos[i] = posList[i].toArray();
                deleChar[i] = charList[i].toArray();
                Arrays.sort(delePos[i]);
            }

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        /**
         *  ################################### taxa 集合 642个taxa
         */
        String vmap2TaxaList = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/taxaList.txt";
        String[] taxa = AoFile.getStringArray(vmap2TaxaList,0);

        double[][] addCount = new double[chrNum][taxa.length];
        int[][] recCount = new int[chrNum][taxa.length];
        int[][] siteWithMinDepthCount = new int[chrNum][taxa.length]; //每个taxa在每条染色体中的有害突变位点

        chrList.parallelStream().forEach(chr -> {
            String delVmapFileS = "chr" + PStringUtils.getNDigitNumber(3, chr) + "_exon_vmap2.1.vcf.gz";
            //开始读写VCF文件
            delVmapFileS = new File(delVCFDirS, delVmapFileS).getAbsolutePath();
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
                    if (temp.startsWith("##")) {
                        continue;
                    }
                    if (temp.startsWith("#CHROM")) {//说明进入taxa信息列
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
                        if (cnt % 1000 == 0) {
                            System.out.println(String.valueOf(cnt) + " lines on chr " + String.valueOf(chr));
                        }
                        List<String> l = PStringUtils.fastSplit(temp.substring(0, 100), "\t");
                        int pos = Integer.valueOf(l.get(1)); //根据pos找到index
                        int index = Arrays.binarySearch(delePos[chrIndex], pos);
                        if (index < 0) {
                            continue; //前面建立的 delePos[][] 和deleChar[][] 都是为现在在vcf文件中找位置贡献的，不是有害突变的位点，都过滤。*************************************************************
                        }
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
                        //code:
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
                            if (depth < minDepth) {
                                continue; //最小测序深度是2，如果小于2，则弃用
                            }
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
                                addCount[chrIndex][taxaIndex] += 0.5;
                            } else {
                                addCount[chrIndex][taxaIndex] += 1;
                                recCount[chrIndex][taxaIndex] += 1;
                            }
                            siteWithMinDepthCount[chrIndex][taxaIndex]++; //taxa有多少个有害突变位点
                        }
                    }
                }
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }

            /**
             * 进行每个样品每条染色体的计算
             */
        });

        try {
            BufferedWriter bw = IOUtils.getTextWriter(addCountFileS);
            bw.write("Taxa\tChr\tDeleteriousCountPerHaplotype\tSiteCountWithMinDepth\tIfSelectedRegion"); //每个taxa有多少个加性效应的derivedAllele 每个taxa在del库中含有基因型的个数
            bw.newLine();
            for (int i = 0; i < addCount.length; i++) {
                String chr = String.valueOf(i + 1);
                for (int j = 0; j < addCount[0].length; j++) {
                    bw.write(taxa[j] + "\t" + chr + "\t" + String.valueOf(addCount[i][j]) + "\t" + String.valueOf(siteWithMinDepthCount[i][j])+"\t1");
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();

            bw = IOUtils.getTextWriter(recCountFileS);
            bw.write("Taxa\tChr\tDeleteriousCountPerHaplotype\tSiteCountWithMinDepth\tSelectedRegion"); //每个taxa有多少个加性效应的derivedAllele 每个taxa在del库中含有基因型的个数
            bw.newLine();
            for (int i = 0; i < addCount.length; i++) {
                String chr = String.valueOf(i + 1);
                for (int j = 0; j < addCount[0].length; j++) {
                    bw.write(taxa[j] + "\t" + chr + "\t" + String.valueOf(recCount[i][j]) + "\t" + String.valueOf(siteWithMinDepthCount[i][j])+"\t1");
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
        }


        /**
         * 将输出的文件添加分组信息：每个taxa的倍性，亚群分布，mutation burden Index
         */

        try{
            String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/taxaList.txt";
            AoFile.readheader(taxaSummaryFileS);
            double depthCut = 3; //只保留深度大于3的taxa样本
            HashMap<String, String> taxaGroupMap = new HashMap();
            HashMap<String, String> taxaSubMap = new HashMap();
            HashMap<String, String> taxaGroupIDMap = new HashMap();
            RowTable t = new RowTable (taxaSummaryFileS);
            ArrayList<String> taxaList = new ArrayList();
            for (int i = 0; i < t.getRowNumber(); i++) {
                taxaGroupMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 10));
                taxaSubMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 11));
                taxaGroupIDMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 7));
                if (t.getCellAsDouble(i, 12) < depthCut) continue;
                taxaList.add(t.getCellAsString(i, 0));
            }
            String[] taxaWithHighDepth = taxaList.toArray(new String[taxaList.size()]); //taxa是具有高深度的taxa列表
            Arrays.sort(taxaWithHighDepth);

            BufferedReader br = AoFile.readFile(addCountFileS);
            BufferedWriter bw = AoFile.writeFile(addCountFileAddGroupS);
            t=new RowTable(addCountFileS);
            List<String> headerl = t.getHeader();
            for (int i = 0; i < headerl.size(); i++) {
                bw.write(headerl.get(i) + "\t");
            }
            bw.write("Group\tSubspecies\tGroupID\tRatio");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) { //这里的t是 addCountFileS
                int index = Arrays.binarySearch(taxa, t.getCellAsString(i, 0));
                if (index < 0) continue;
                if(taxaGroupMap.get(taxa[index]).equals("ExclusionHexaploid") || taxaGroupMap.get(taxa[index]).equals("ExclusionTetraploid")) continue;
                double genotypesite = Double.valueOf(t.getCellAsDouble(i, 2));
                if(genotypesite == 0) continue;
                double ratio = Double.valueOf(t.getCellAsDouble(i, 2))/Double.valueOf(t.getCellAsDouble(i, 3));
                for (int j = 0; j < t.getColumnNumber(); j++) { //按列书写
                    bw.write(t.getCellAsString(i,j) + "\t");
                }
                bw.write(taxaGroupMap.get(taxa[index])
                        +"\t"+taxaSubMap.get(taxa[index])
                        +"\t"+taxaGroupIDMap.get(taxa[index])
                        +"\t"+String.format("%.4f",ratio));
                bw.newLine();
            }
            bw.close();
        }
        catch(Exception e){

        }
    }


    public void mergeExonSNPAnnotation(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/003_exonSNPAnnotation";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/004_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
        AoFile.mergeTxt(infileDirS,outfileS);

    }

}
