package WheatGeneticLoad;

import AoUtils.AoFile;
import AoUtils.Bin;
import AoUtils.CalVCF;
import AoUtils.CountSites;
import gnu.trove.list.array.TCharArrayList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.utils.wheat.RefV1Utils;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.table.RowTable;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

public class EstSFS {

    public EstSFS(){
//        this.sampleOutdata();
//        this.addAnc();

//        this.mergeDelgenicSNPAnnotation();
//        this.countDeleteriousVMapII_byChr();
//        this.countDeleteriousVMapIIHighDepth();
//        this.mergeFinalfilebySub();

//        this.getNonsynonymousGerp1();
//        this.diffcheck("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/003_ancestralAllele", "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/003_ancestral/hv_brdis/002_byChrID");


//        this.checkAncestraldiff();

        /**
         * 处理康李鹏的ancestral allele barley secale (maximum likelihood method)
         */
//        this.splitChrfromAncestralLipeng();
//        this.addAnc();
//        this.mergeExonSNPAnnotation();

        /**
         * 计算DAF值，并且生成表格用来画图
         */
//        this.calDAF(); //单线程计算
//        this.addDAF_parallel(); //多线程运行
//        this.mergeExonSNPAnnotation(); //计算后结果合并
//        this.mkDAFtable(); //开始分bin

        /**
         * 处理达兴的简约法ancestral allele Pasimony
         */
//        this.extractAncestralAllele();
//        this.test();
//        this.deleteriousCount();

//        this.getjDAFtable();
//        this.getjDAFtable_fromrealPop();

        /**
         * 判断过滤heter后的maf分布
         */

//        this.getMAFfromAetauschii();
//        this.mergeTxt();
//        this.maf();
//        CountSites.countSites_singleStream("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/002_exonSNPVCF");

//        CountSites.countSites_singleStream("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/010_exonSNPVCF_filterHeter0.05");

    }

    public void maf(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/007_ancestral_Barley_secale_parsimony/007_maf/002_merge001";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/007_ancestral_Barley_secale_parsimony/007_maf/003_MAFtable";
        int bins = 20;
        int columIndex = 0;

        Bin.mkBarplotofMAF(infileDirS,outfileDirS,bins,columIndex);
    }

    public void mergeTxt(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/007_ancestral_Barley_secale_parsimony/007_maf/001_DDmaf";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/007_ancestral_Barley_secale_parsimony/007_maf/002_merge001/ae.tauschii_exon_vmap2.1_filterbyHeter0.05_maf.txt.gz";
        new CountSites().mergeTxt(infileDirS,outfileS);
    }

    public void getMAFfromAetauschii(){ //本地运行常用
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/010_exonSNPVCF_filterHeter0.05";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/007_ancestral_Barley_secale_parsimony/007_maf/001_DDmaf";
        String[] chrArr ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};
        String taxaListS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/002_groupbyPloidy_removeBadTaxa/Ae.tauschii_S36.txt";


        for (int j = 0; j < chrArr.length; j++) {
            String infileS = new File(infileDirS,"chr" + chrArr[j] + "_exon_vmap2.1_filterbyHeter0.05.vcf.gz").getAbsolutePath();
            String outfileS = new File(outfileDirS,"chr" + chrArr[j] + "_exon_vmap2.1_filterbyHeter0.05_maf.txt.gz").getAbsolutePath();
            CalVCF.calMafFromPop(infileS,outfileS,taxaListS);
        }
    }

    /**
     * 真正计算六倍体四倍体的群体数目，建立表格
     * 六倍体：419
     * 四倍体：187
     */
    public void getjDAFtable_fromrealPop(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/009_exonSNPAnnotation_addAnc_addDAF_barley_secalePasimony/chr001_SNP_anno.txt.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/006_ancestralfromLipeng/004_DAFtable/jDAF/jDAFchr001_pop.txt";
        TDoubleArrayList dafList1 = new TDoubleArrayList();
        TDoubleArrayList dafList2 = new TDoubleArrayList();
        int bins = 20;
        int bin1 = 419*2+1;
        int bin2 = 187*2+1;
        int[][] out = new int[bin1][bin2];
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            AoFile.readheader(infileS);
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                int chrID = Integer.parseInt(l.get(1));
                int posID = Integer.parseInt(l.get(2));
                String chr = RefV1Utils.getChromosome(chrID,posID);
//                    if(chr.contains("D"))continue; //只能用于AABB总体画图
//                    if(chr.contains("A") || chr.contains("B"))continue; //只能用于DD总体画图

//                    if(chr.contains("B") || chr.contains("D"))continue; //只能用于A亚基因组
//                    if(chr.contains("A") || chr.contains("D"))continue; //只能用于B亚基因组
//                if(chr.contains("A") || chr.contains("B"))continue; //只能用于D亚基因组
                System.out.println(temp);
                String type = l.get(12);
                String siftscore = l.get(13);
                String gerpscore = l.get(20);
//                    String phylopscore = l.get(18);
                String daf1 = l.get(33);
                String daf2 = l.get(34);
                if (daf1.startsWith("N") || daf2.startsWith("N")) continue;
                dafList1.add(Double.parseDouble(daf1));
                dafList2.add(Double.parseDouble(daf2));
            }
            br.close();

            double length = Double.parseDouble("1");
            //先建立bound数组,假如 bin=20,那么每个bin的大小是0.05 bound[0] =0; bound[1]=0.05, bound[19]=0.95 即只取了区间的左边开始位置
            double[] bound = new double[bins];
            for (int i = 1; i < bound.length; i++) {
                bound[i] = (double) length / bins * i;
            }

            double[] daf = new double[bins]; //确定落入每个区间的个数，进行计算，20个bin，20个count数字
            for (int i = 0; i < dafList1.size(); i++) { //对每个daf值进行区间的判断
                double value = dafList1.get(i);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) {
                    index = -index - 2;
                }
                for (int j = 0; j < dafList2.size(); j++) {
                    double value2 = dafList2.get(j);
                    int index2 = Arrays.binarySearch(bound, value2);
                    if (index2 < 0) {
                        index2 = -index2 - 2;
                    }
                    out[index][index2]++;
                }
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
            }

            for (int i = 0; i < bound.length; i++) {
                String coordinate = String.format("%.3f", (double) bound[i] + (double) (length / bins) / (double) 2);
                StringBuilder sb = new StringBuilder();
                if (i==0){
                    bw.write(coordinate);
                    continue;
                }
                sb.append("\t").append(coordinate);
                bw.write(sb.toString());
            }
            bw.newLine();

            for (int i = 0; i < out.length; i++) {
                String coordinate = String.format("%.3f", (double) bound[i] + (double) (length / bins) / (double) 2);
                bw.write(coordinate);
                for (int j = 0; j < out[0].length; j++) {
                    bw.write("\t" + out[i][j]);
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }



    /**
     *
     */
    public void getjDAFtable(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/009_exonSNPAnnotation_addAnc_addDAF_barley_secalePasimony/chr001_SNP_anno.txt.gz";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/009_exonSNPAnnotation_addAnc_addDAF_barley_secalePasimony/chr005_SNP_anno.txt.gz";

        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/006_ancestralfromLipeng/004_DAFtable/jDAF/jDAFchr001.txt";
        TDoubleArrayList dafList1 = new TDoubleArrayList();
        TDoubleArrayList dafList2 = new TDoubleArrayList();
        int bins = 20;
        int[][] out = new int[20][20];
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            AoFile.readheader(infileS);
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                int chrID = Integer.parseInt(l.get(1));
                int posID = Integer.parseInt(l.get(2));
                String chr = RefV1Utils.getChromosome(chrID,posID);
//                    if(chr.contains("D"))continue; //只能用于AABB总体画图
//                    if(chr.contains("A") || chr.contains("B"))continue; //只能用于DD总体画图

//                    if(chr.contains("B") || chr.contains("D"))continue; //只能用于A亚基因组
//                    if(chr.contains("A") || chr.contains("D"))continue; //只能用于B亚基因组
//                if(chr.contains("A") || chr.contains("B"))continue; //只能用于D亚基因组
                System.out.println(temp);
                String type = l.get(12);
                String siftscore = l.get(13);
                String gerpscore = l.get(20);
//                    String phylopscore = l.get(18);
                String daf1 = l.get(33);
                String daf2 = l.get(34);
//
//                String daf1 = l.get(8);
//                String daf2 = l.get(9);
                if (daf1.startsWith("N") || daf2.startsWith("N")) continue;
                dafList1.add(Double.parseDouble(daf1));
                dafList2.add(Double.parseDouble(daf2));
            }
            br.close();

            double length = Double.parseDouble("1");
            //先建立bound数组,假如 bin=20,那么每个bin的大小是0.05 bound[0] =0; bound[1]=0.05, bound[19]=0.95 即只取了区间的左边开始位置
            double[] bound = new double[bins];
            for (int i = 1; i < bound.length; i++) {
                bound[i] = (double) length / bins * i;
            }

            double[] daf = new double[bins]; //确定落入每个区间的个数，进行计算，20个bin，20个count数字
            for (int i = 0; i < dafList1.size(); i++) { //对每个daf值进行区间的判断
                double value = dafList1.get(i);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) {
                    index = -index - 2;
                }
                for (int j = 0; j < dafList2.size(); j++) {
                    double value2 = dafList2.get(j);
                    int index2 = Arrays.binarySearch(bound, value2);
                    if (index2 < 0) {
                        index2 = -index2 - 2;
                    }
                    out[index][index2]++;
                }
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
            }

            for (int i = 0; i < bound.length; i++) {
//                String coordinate = String.format("%.3f", (double) bound[i] + (double) (length / bins) / (double) 2);
                String coordinate = String.format("%.3f", (double) bound[i]*20);
                StringBuilder sb = new StringBuilder();
                if (i==0){
                    bw.write(coordinate);
                    continue;
                }
                sb.append("\t").append(coordinate);
                bw.write(sb.toString());
            }
            bw.newLine();

            for (int i = 0; i < out.length; i++) {
//                String coordinate = String.format("%.3f", (double) bound[i] + (double) (length / bins) / (double) 2);
                String coordinate = String.format("%.3f", (double) bound[i]*20);

                bw.write(coordinate);
                for (int j = 0; j < out[0].length; j++) {
                    bw.write("\t" + out[i][j]);
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void deleteriousCount(){
//        new DeleteriousCountbyPop().countDeleteriousVMapII_byChr();
//        new DeleteriousCountbyPop().DeltoSynonymousRatio();
    }

    public void test (){
        String SNPAnnoFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/004_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz"; //有害变异信息库
        AoFile.readheader(SNPAnnoFileS);
        System.out.println("begin");

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
                int index = Integer.parseInt(l.get(1)) - 1; //染色体号的索引
                int pos = Integer.parseInt(l.get(2));
                String variantType = l.get(12);
                String sift = l.get(13);
                String gerp = l.get(20);
                /**
                 * 定义有害突变，不是有害突变，就忽略不计 ################ 需要修改 需要修改 需要修改 ################
                 */

//                if (!variantType.equals("NONSYNONYMOUS") || sift.equals("NA") || gerp.equals("NA"))continue;
//                double siftd = Double.parseDouble(sift);
//                double gerpd = Double.parseDouble(gerp);
//                if (siftd >= 0.05 || gerpd <= 1)continue;

                /**
                 * 定义有害突变，不是有害突变，就忽略不计 ################ 需要修改 需要修改 需要修改 ################
                 */
                if (!variantType.equals("SYNONYMOUS"))continue;

                //################### 需要修改 //###################//###################//###################//###################
//                String ancestralAllele = l.get(22); //不同的数据库，这一列的信息不一样，千万要注意!!!!!!!!!!!!!!!!! 祖先基因的数据库
//                String ancestralAllele = l.get(15); //不同的数据库，这一列的信息不一样，千万要注意!!!!!!!!!!!!!!!!! 祖先基因的数据库
                String ancestralAllele = l.get(31);
                //################### 需要修改 //###################//###################//###################//###################

                String majorAllele = l.get(5);
                String minorAllele = l.get(6);
                if (ancestralAllele.equals(majorAllele)) {
                    derivedAllele = minorAllele;
                    posList[index].add(Integer.parseInt(l.get(2))); //将包含有derived allele的位点添加到Poslist
                    charList[index].add(derivedAllele.charAt(0)); //返回的是char类型的字符
                }
                if (ancestralAllele.equals(minorAllele)) {
                    derivedAllele = majorAllele;
                    posList[index].add(Integer.parseInt(l.get(2)));
                    charList[index].add(derivedAllele.charAt(0));
                }
                else if (!(ancestralAllele.equals(majorAllele) || ancestralAllele.equals(minorAllele))){
                }
                cnt++;
                if (cnt%1000==0) System.out.println("cnt is " + cnt);

            }
            br.close();
            System.out.println();

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
    }

    /**
     * 提取达兴的简约法的ancestral allele
     */
    public void extractAncestralAllele(){
        AoFile.readheader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/007_exonSNPAnnotation_addDAF_barleyUratu/chr001_SNP_anno.txt.gz");
        System.out.println("*****************************");
        AoFile.readheader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/008_exonSNPAnnotation_addAnc/chr001_SNP_anno.txt.gz");
        String infileDirS ="/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/007_exonSNPAnnotation_addDAF_barleyUratu";
        String dbDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/007_ancestral_Barley_secale_parsimony/002_byChrID";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/008_exonSNPAnnotation_addAnc";
        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
//        fsList.stream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String dbS = f.getName().substring(0,6) + "_secer_hv_ancestral.txt.gz";
                dbS = new File(dbDirS,dbS).getAbsolutePath();
                String outfileS = new File(outfileDirS,f.getName()).getAbsolutePath();
                BufferedReader br = AoFile.readFile(dbS);
                String temp = null;
                String header = br.readLine();
                List<String> l = new ArrayList<>();
                TIntArrayList posl = new TIntArrayList();
                TCharArrayList ancl = new TCharArrayList();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(1));
                    String anc = l.get(2);
                    posl.add(pos);
                    ancl.add(anc.charAt(0));
                }
                br.close();
                System.out.println(f.getName() + "\t DB is completed.");
                int[] posArray = posl.toArray(new int[posl.size()]);
                char[] ancArray = ancl.toArray(new char[ancl.size()]);
                Arrays.sort(posArray);


                br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                header = br.readLine();
                bw.write(header + "\tAncestral_barley_secale_parsimony");
                bw.newLine();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    int pos = Integer.parseInt(l.get(2));
                    int index = Arrays.binarySearch(posArray,pos);
                    if (index > -1){
                        bw.write(temp + "\t" + ancArray[index]);
                        bw.newLine();
                    }else if (index <0){
                        bw.write(temp + "\tNA");
                        bw.newLine();
                    }
                }
                br.close();
                bw.flush();
                bw.close();
                System.out.println(f.getName() + "\tis completed.");

            } catch (Exception e) {
                e.printStackTrace();
            }
        });

    }

    public void mkDAFtable(){
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/004_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
//        AoFile.readheader(infileS);
        new Bin().getDAFtable();

//

    }

    public void runParallele_listFile(){ //本地运行常用
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/005_exonSNPAnnotation";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/006_exonSNPAnnotation_addDAF";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/006_exonSNPAnnotation_addDAF";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/007_exonSNPAnnotation_addDAF_barleyUratu";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/008_exonSNPAnnotation_addAnc";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/009_exonSNPAnnotation_addAnc_addDAF_barley_secalePasimony";
        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        fsList.stream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = new File(outfileDirS,f.getName()).getAbsolutePath();
            this.calDAF(f.getAbsolutePath(),outfileS);
            System.out.println(f.getName() + "\tis completed at " + outfileS);
        });
    }

    /**
     * Goal:根据康李鹏提供的 ancestral allele，计算Daf,Daf_ABD Daf_AB Daf_D
     */
    public void calDAF(String dbfileS, String outfileS) { //String dbfileS, String ancS, String outfileS
//        String dbfileS = "";
//        String outfileS = "";
        boolean ifd = false;
        double daf = Double.NaN;
        double daf_ABD = Double.NaN;
        double daf_AB = Double.NaN;
        int cntAncNum = 0;
        File f = new File(dbfileS); //根据ancestral allele 文件，得到染色体号
        int chr = Integer.parseInt(f.getName().substring(3, 6));
        //根据染色体号进行AB还是D的判断
        int[] db = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
        Arrays.sort(db);
        if (Arrays.binarySearch(db, chr) > -1) { //说明是属于D的
            ifd = true;
        }
        try {
            String chrS = PStringUtils.getNDigitNumber(3, chr);
            BufferedReader br = AoFile.readFile(dbfileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = br.readLine(); //read header
            if (ifd == false) {
//                bw.write(temp + "\tDaf_barley_secale\tDaf_ABD_barley_secale\tDaf_AB_barley_secale");
//                bw.write(temp + "\tDaf_barley_urartu\tDaf_ABD_barley_urartu\tDaf_AB_barley_urartu");
                bw.write(temp + "\tDaf_barley_secalePasimony\tDaf_ABD_barley_secalePasimony\tDaf_AB_barley_secalePasimony");
                bw.newLine();
            } else if (ifd == true) {
//                bw.write(temp + "\tDaf_barley_secale\tDaf_ABD_barley_secale\tDaf_D_barley_secale");
//                bw.write(temp + "\tDaf_barley_urartu\tDaf_ABD_barley_urartu\tDaf_AB_barley_urartu");
                bw.write(temp + "\tDaf_barley_secalePasimony\tDaf_ABD_barley_secalePasimony\tDaf_AB_barley_secalePasimony");
                bw.newLine();
            }

            int cntAnc = 0;
            int cntAncNotMajororMinor = 0;
            while ((temp = br.readLine()) != null) {
                List<String> l = PStringUtils.fastSplit(temp);
                int pos = Integer.parseInt(l.get(2));
                String major = l.get(5);
                String minor = l.get(6);
                double maf = Double.parseDouble(l.get(7));
                double AAF_ABD = Double.parseDouble(l.get(8));
                double AAF_AB = Double.parseDouble(l.get(9));

                //################### 需要修改 //###################//###################//###################//###################
//               String ancAllele = l.get(22); //不同的数据库，这一列的信息不一样，千万要注意!!!!!!!!!!!!!!!!! 祖先基因的数据库
//                String ancAllele = l.get(15); //不同的数据库，这一列的信息不一样，千万要注意!!!!!!!!!!!!!!!!! 祖先基因的数据库
               String ancAllele = l.get(31);
                //################### 需要修改 //###################//###################//###################//###################

                StringBuilder sb = new StringBuilder();
                if (!ancAllele.equals("NA")) { //表明含有anc
                    //如果ancestral allele存在,且等于major，则derived allele等于minor, daf 就等于maf
                    //如果ancestral allele存在,且等于minor，则derived allele等于major, daf 就等于 1-daf1
                    if (ancAllele.equals(minor)) {
                        cntAnc++;
                        daf = 1 - maf;
                        if (AAF_ABD > 0.5) { //说明AAF_ABD是major， 祖先状态是minor的，所有DAF是major
                            daf_ABD = AAF_ABD;
                        } else if (AAF_ABD < 0.5) { //说明AAF_ABD是minor， 祖先状态是minor的，所有DAF是major
                            daf_ABD = 1 - AAF_ABD;
                        }
                        if (AAF_AB > 0.5) {
                            daf_AB = AAF_AB;
                        } else if (AAF_AB < 0.5) {
                            daf_AB = 1 - AAF_AB;
                        }
                        //多加一道判断，如果群体内部没有分离，直接将DAF设置为NA
                        String DAF_ABD = String.format("%.4f", daf_ABD);
                        String DAF_AB = String.format("%.4f", daf_AB);
                        if (DAF_ABD.equals("0.0000") || DAF_ABD.equals("1.0000")) {
                            DAF_ABD = "NA";
                        }
                        if (DAF_AB.equals("0.0000") || DAF_AB.equals("1.0000")) {
                            DAF_AB = "NA";
                        }
                        sb.append(temp).append("\t").append(String.format("%.4f", daf)).append("\t").append(DAF_ABD).append("\t").append(DAF_AB);
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    if (ancAllele.equals(major)) {
                        cntAnc++;
                        daf = maf;
                        if (AAF_ABD > 0.5) { //说明AAF_ABD是major， 祖先状态是major的，所有DAF是minor
                            daf_ABD = 1 - AAF_ABD;
                        } else if (AAF_ABD < 0.5) {
                            daf_ABD = AAF_ABD;
                        }
                        if (AAF_AB > 0.5) {
                            daf_AB = 1 - AAF_AB;
                        } else if (AAF_AB < 0.5) {
                            daf_AB = AAF_AB;
                        }
                        //多加一道判断，如果群体内部没有分离，直接将DAF设置为NA
                        String DAF_ABD = String.format("%.4f", daf_ABD);
                        String DAF_AB = String.format("%.4f", daf_AB);
                        if (DAF_ABD.equals("0.0000") || DAF_ABD.equals("1.0000")) {
                            DAF_ABD = "NA";
                        }
                        if (DAF_AB.equals("0.0000") || DAF_AB.equals("1.0000")) {
                            DAF_AB = "NA";
                        }
                        sb.append(temp).append("\t").append(String.format("%.4f", daf)).append("\t").append(DAF_ABD).append("\t").append(DAF_AB);
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    if (!ancAllele.equals(minor) && (!ancAllele.equals(major))) {
                        sb.append(temp).append("\t").append("NA").append("\t").append("NA").append("\t").append("NA");
                        bw.write(sb.toString());
                        bw.newLine();
                        //System.out.println("CHR" + PStringUtils.getNDigitNumber(3, CHR) + "\t" + pos + " are neither major nor minor.");
                        cntAncNotMajororMinor++;
                    }

                } else { //表明不含anc
                    sb.append(temp).append("\t").append("NA").append("\t").append("NA").append("\t").append("NA");
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            double ratio = (double) cntAncNotMajororMinor / (cntAncNotMajororMinor + cntAnc);
            bw.flush();
            bw.close();
            br.close();
            System.out.println(f.getName() + "\tis completed at " + outfileS + "\t" + cntAnc + "\tancestral allele are with daf value by state major or minor");
            System.out.println(new File(dbfileS).getName() + "\thave " + cntAncNotMajororMinor + " sites which are neither major nor minor. The ratio is " + String.format("%.4f", ratio));
        } catch (Exception e) {
            e.printStackTrace();
        }
    }




    public void mergeExonSNPAnnotation(){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/007_exonSNPAnnotation_addDAF_barleyUratu";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/008_exonSNPAnnotation_addAnc";
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/009_exonSNPAnnotation_addAnc_addDAF_barley_secalePasimony";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/004_exonSNPAnnotation_merge/001_exonSNP_anno.txt.gz";
        AoFile.mergeTxt(infileDirS,outfileS);
    }

    /**
     * 将康李鹏的结果进行拆分,添加到SNP数据库中
     */
    public void splitChrfromAncestralLipeng(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/006_ancestralfromLipeng/001_ori/wheat.anc_1.gz";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/006_ancestralfromLipeng/002_byChrID";
        int chrNum = 42;
        try {
            String[] outfileS = new String[chrNum];
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter[] bw = new BufferedWriter[chrNum];
            //确定每个文件的输出路径，每个文件都写一个表头chr pos ancestral
            for (int i = 0; i < chrNum; i++) {
                String chr = PStringUtils.getNDigitNumber(3,i+1);
                outfileS[i] = new File(outfileDirS,"chr" + chr + "_barleyVSsecale_ancestralAllele.txt.gz").getAbsolutePath();
                bw[i] = AoFile.writeFile(outfileS[i]);
                bw[i].write("Chr\tPos\tAncestralAllele");
                bw[i].newLine();
            }

            String temp = null;
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#"))continue;
                l = PStringUtils.fastSplit(temp);
                String chrom = l.get(0).substring(3);
                int chrpos = Integer.parseInt(l.get(2));
                String anc = l.get(10);
                int chr = RefV1Utils.getChrID(chrom,chrpos);
                int pos = RefV1Utils.getPosOnChrID(chrom,chrpos);
                int index = chr-1;
                bw[index].write(chr + "\t" + pos+"\t"+anc);
                bw[index].newLine();
            }
            br.close();
            for (int i = 0; i < chrNum; i++) {
                bw[i].flush();
                bw[i].close();
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     * 为了查看 2个outgroup和3个outgroup计算出来的ancestral allele 之间的区别
     */
    public void checkAncestraldiff(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/000_/del";
        String outfileDirS = "";
        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                new AoFile().readheader(infileS);
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "_subset.txt.gz").getAbsolutePath();
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().split(".txt.gz")[0] + "_subset.txt.gz").getAbsolutePath();
                }
                String temp = null;
                List<String> l = new ArrayList<>();
                int cnt=0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    cnt++;
                    String chr = l.get(1);
                    String an1 = l.get(22);
                    String an2 = l.get(23);
                    if(!an1.equals(an2)){
                        System.out.println(an1 + "\t" + an2 + "\t" + cnt + "\t" + chr);
                    }

                }
                br.close();
                System.out.println(f.getName() + "\tis completed at ");
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    /**
     * 查看达兴2次给的ancestral allele （2个外围群的值）md5为什么不一样
     */
    public void diffcheck(){
        String inputDir1 = "";
        String inputDir2 = "";
        List<File> files1=IOUtils.getVisibleFileListInDir(inputDir1);
        List<File> files2=IOUtils.getVisibleFileListInDir(inputDir2);
        BufferedReader br;
        BufferedReader br2;
        for (int i = 0; i < files1.size(); i++) {
            try {
                br = IOUtils.getTextGzipReader(files1.get(i).getAbsolutePath());
                br2 = IOUtils.getTextGzipReader(files2.get(i).getAbsolutePath());
                String temp = null;
                List<String> l = new ArrayList<>();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    String temp2 = br2.readLine();
                    if(!temp.equals(temp2)){
                        System.out.println(temp + "\t" + temp2);
                    }

                }
                br.close();
                br2.close();
                System.out.println();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        }



    }

    public void getNonsynonymousGerp1(){
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/002_/del_merged/chr_SNP_anno_del.txt.gz";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/005_nonsynonymousGerp1/del_merged";

        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/003_/syn_merge/chr_SNP_anno_syn.txt.gz";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/006_synonymousGerp1/syn_merge";


        new AoFile().readheader(infileS);
        try {
            String outfileS = new File(outfileDirS,new File(infileS).getName()).getAbsolutePath();
            BufferedReader br = IOUtils.getTextGzipReader(infileS);
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            String temp = br.readLine();
            bw.write(temp);
            bw.newLine();
            List<String> l;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                if(l.get(20).startsWith("N"))continue;
                double gerp = Double.parseDouble(l.get(20));
                if(gerp>1){
                    bw.write(temp);
                    bw.newLine();
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



    //根据最终生成的文件，进行 A B D sub的合并
    public void mergeFinalfilebySub(){
        //change
//        String infileS = "";

//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/001_/additiveDeleterious_vmap2_highDepth_bychr.txt";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/002_/additiveDeleterious_vmap2_highDepth_bychr.txt";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/003_/additiveDeleterious_vmap2_highDepth_bychr.txt";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/000_/additiveDeleterious_vmap2_highDepth_bychr.txt";


//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/002_delCount_highDepth/000_/additiveDeleterious_vmap2_highDepth_bychr.txt";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/002_delCount_highDepth/001_/additiveDeleterious_vmap2_highDepth_bychr.txt";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/002_delCount_highDepth/002_/additiveDeleterious_vmap2_highDepth_bychr.txt";
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/002_delCount_highDepth/003_/additiveDeleterious_vmap2_highDepth_bychr.txt";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/002_delCount_highDepth/005_nonsynonymousGerp1/additiveDeleterious_vmap2_highDepth_bychr.txt";
//        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/002_delCount_highDepth/006_synonymousGerp1/additiveDeleterious_vmap2_highDepth_bychr.txt";


        String splitDirS = new File(infileS).getParent()+"/split";
        new File(splitDirS).mkdirs();

        String outfileS = infileS.replaceFirst("_bychr.txt", "_bysub.txt");

        //查看有多少个taxa
        RowTable<String> t = new RowTable<>(infileS);
        Set<String> s = new HashSet<>(t.getColumn(0));
        String[] taxa = s.toArray(new String[s.size()]);
        Arrays.sort(taxa);
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String header = br.readLine();
            BufferedWriter[] bw = new BufferedWriter[taxa.length];
            for (int i = 0; i < taxa.length; i++) {
                String tempoutS = new File(splitDirS,taxa[i]+".txt").getAbsolutePath();
                bw[i] = IOUtils.getTextWriter(tempoutS);
                bw[i].write(header);
                bw[i].newLine();
            }

            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String ta = l.get(0);
                int index = Arrays.binarySearch(taxa,ta);
                bw[index].write(temp);
                bw[index].newLine();
                cnt++;

            }
            br.close();
            for (int i = 0; i < taxa.length; i++) {
                bw[i].flush();
                bw[i].close();
            }
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }


        //不变的
        String taxaGroupFileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/001_toFeiLu/004_wheatVMapII_germplasmInfoAddSubspecies_20191225.txt";
        String taxaGroup2FileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";
        t = new RowTable (taxaGroupFileS);
        HashMap<String, String> taxaGroupMap = new HashMap();
        HashMap<String, String> taxaSubMap = new HashMap();
        HashMap<String, String> taxaGroupIDMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaGroupMap.put(t.getCellAsString(i, 4), t.getCellAsString(i, 26));
            taxaSubMap.put(t.getCellAsString(i, 4), t.getCellAsString(i, 27));
        }
        t = new RowTable (taxaGroup2FileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaGroupIDMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 7));
        }

        HashMap<String,Integer> h = new HashMap<>();
        h.put("A",0);
        h.put("B",1);
        h.put("D",2);

        String outDirS = null;
        for (int i = 0; i < taxa.length; i++) {
            try {
                outDirS = new File(splitDirS).getParent() + "/split_merge";
                new File(outDirS).mkdirs();

                String inS = new File(splitDirS,taxa[i]+".txt").getAbsolutePath();
                String outS = new File(outDirS,taxa[i]+"_bysub.txt").getAbsolutePath();
                BufferedReader br = IOUtils.getTextReader(inS);
                BufferedWriter bw = IOUtils.getTextWriter(outS);
                bw.write(br.readLine().replaceFirst("Chr","Sub"));
                bw.newLine();
                TDoubleArrayList[] derivedDelList = new TDoubleArrayList[3];
                TDoubleArrayList[] genotypeList = new TDoubleArrayList[3];
                for (int j = 0; j < derivedDelList.length; j++) {
                    derivedDelList[j] = new TDoubleArrayList();
                    genotypeList[j] = new TDoubleArrayList();

                }
                String temp = null;
                List<String> l = new ArrayList<>();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    cnt++;
                    int chr = Integer.parseInt(l.get(1));
                    double d = Double.parseDouble(l.get(2));
                    double k = Double.parseDouble(l.get(3));
                    String sub = this.convertChrtoSub(chr);
                    derivedDelList[h.get(sub)].add(d);
                    genotypeList[h.get(sub)].add(k);
                }
                br.close();


                //求和
                DescriptiveStatistics[] d = new DescriptiveStatistics[3];
                DescriptiveStatistics[] dd = new DescriptiveStatistics[3];

                Double[] ratio = new Double[3];
                for (int j = 0; j < derivedDelList.length; j++) {
                    d[j] = new DescriptiveStatistics(derivedDelList[j].toArray());
                    dd[j] = new DescriptiveStatistics(genotypeList[j].toArray());
                    ratio[j] = d[j].getSum()/dd[j].getSum();

                }

                HashMap<Integer,String> hhh = new HashMap<>();
                hhh.put(0,"A");
                hhh.put(1,"B");
                hhh.put(2,"D");


                for (int j = 0; j < derivedDelList.length; j++) {
                    if(dd[j].getSum()==0)continue;
                    bw.write(taxa[i] + "\t" + hhh.get(j) + "\t" + String.format("%.1f",d[j].getSum()) + "\t" + String.format("%.0f",dd[j].getSum())
                            + "\t" + taxaGroupMap.get(taxa[i]) + "\t" + taxaSubMap.get(taxa[i]) + "\t" + taxaGroupIDMap.get(taxa[i]) + "\t" + String.format("%.3f",ratio[j]));
                    bw.newLine();
                }

                bw.flush();
                bw.close();

            } catch (Exception e) {
                e.printStackTrace();
            }

        }

        new CountSites().mergeTxtbysuffix(outDirS,outfileS,".txt");
        new File(splitDirS).delete();
        new File(outDirS).delete();

    }

    /**
     * 根据染色体号，返回所在亚基因组
     * @param a
     * @return
     */
    public String convertChrtoSub(int a){
        int[] arra = {1, 2, 7, 8, 13, 14, 19, 20, 25, 26, 31, 32, 37, 38};
        int[] arrb = {3, 4, 9, 10, 15, 16, 21, 22, 27, 28, 33, 34, 39, 40};
        int[] arrd = {5, 6, 11, 12, 17, 18, 23, 24, 29, 30, 35, 36, 41, 42};
        HashMap<Integer, String> hml = new HashMap<>();
        Arrays.sort(arra);
        Arrays.sort(arrb);
        Arrays.sort(arrd);
        for (int i = 0; i < arra.length; i++) {
            hml.put(arra[i], "A");
            hml.put(arrb[i], "B");
            hml.put(arrd[i], "D");
        }
        String subgenome = hml.get(a);
        return subgenome;
    }



    private void countDeleteriousVMapIIHighDepth () {
        //不变的文件
        String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxaDepth_summary.txt";
        String taxaGroupFileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/001_toFeiLu/004_wheatVMapII_germplasmInfoAddSubspecies_20191225.txt";
        String taxaGroup2FileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";

        //可变的文件
//        String addInfileS = "";
//        String addOutfileS = "";

//        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/002_/additiveDeleterious_vmap2_bychr.txt";
//        String addOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/002_/additiveDeleterious_vmap2_highDepth_bychr.txt";

//        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/003_/additiveDeleterious_vmap2_bychr.txt";
//        String addOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/003_/additiveDeleterious_vmap2_highDepth_bychr.txt";

//        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/000_/additiveDeleterious_vmap2_bychr.txt";
//        String addOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/000_/additiveDeleterious_vmap2_highDepth_bychr.txt";



//        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/001_delCount/000_/additiveDeleterious_vmap2_bychr.txt";
//        String addOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/002_delCount_highDepth/000_/additiveDeleterious_vmap2_highDepth_bychr.txt";

//        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/001_delCount/001_/additiveDeleterious_vmap2_bychr.txt";
//        String addOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/002_delCount_highDepth/001_/additiveDeleterious_vmap2_highDepth_bychr.txt";

//        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/001_delCount/002_/additiveDeleterious_vmap2_bychr.txt";
//        String addOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/002_delCount_highDepth/002_/additiveDeleterious_vmap2_highDepth_bychr.txt";

        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/001_delCount/003_/additiveDeleterious_vmap2_bychr.txt";
        String addOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/002_delCount_highDepth/003_/additiveDeleterious_vmap2_highDepth_bychr.txt";

//        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/001_delCount/004_gerp/additiveDeleterious_vmap2_bychr.txt";
//        String addOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/002_delCount_highDepth/004_gerp/additiveDeleterious_vmap2_highDepth_bychr.txt";

//        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/001_delCount/005_nonsynonymousGerp1/additiveDeleterious_vmap2_bychr.txt";
//        String addOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/002_delCount_highDepth/005_nonsynonymousGerp1/additiveDeleterious_vmap2_highDepth_bychr.txt";

//        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/001_delCount/006_synonymousGerp1/additiveDeleterious_vmap2_bychr.txt";
//        String addOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/002_delCount_highDepth/006_synonymousGerp1/additiveDeleterious_vmap2_highDepth_bychr.txt";





//        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/additiveDeleterious_vmap2.txt";
//        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/001_/additiveDeleterious_vmap2_bychr.txt";
//        String addOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/001_/additiveDeleterious_vmap2_highDepth_bychr.txt";
//        String recInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/001_/recessiveDeleterious_vmap2_bychr.txt";
//        String recOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/001_/recessiveDeleterious_vmap2_highDepth_bychr.txt";
//        String taxaGroupFileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/001_toFeiLu/004_wheatVMapII_germplasmInfoAddSubspecies_20191225.txt";
//        String taxaGroup2FileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";


//        String taxaSummaryFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/002_VMapIIDepth/taxaDepth_summary.txt";
//        String addInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/additiveSyn_vmap2.txt";
//        String addOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/additiveSyn_vmap2_highDepth.txt";
//        String recInfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/recessiveSyn_vmap2.txt";
//        String recOutfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/003_delCount_highDepth/recessiveSyn_vmap2_highDepth.txt";
//        String taxaGroupFileS = "/Users/Aoyue/project/wheatVMapII/001_germplasm/GermplasmDB/001_toFeiLu/004_wheatVMapII_germplasmInfoAddSubspecies_20191225.txt";
//        String taxaGroup2FileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/006_tree/005_ABsub_maf0.01_20191207/000_prepareData/001_input/taxaList.txt";



        new AoFile().readheader(taxaGroupFileS);
        double depthCut = 3;
        RowTable t = new RowTable (taxaSummaryFileS);
        ArrayList<String> taxaList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.getCellAsDouble(i, 2) < depthCut) continue;
            taxaList.add(t.getCellAsString(i, 0));
        }
        String[] taxa = taxaList.toArray(new String[taxaList.size()]); //taxa是具有高深度的taxa列表
        Arrays.sort(taxa);
        t = new RowTable (taxaGroupFileS);
        HashMap<String, String> taxaGroupMap = new HashMap();
        HashMap<String, String> taxaSubMap = new HashMap();
        HashMap<String, String> taxaGroupIDMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaGroupMap.put(t.getCellAsString(i, 4), t.getCellAsString(i, 26));
            taxaSubMap.put(t.getCellAsString(i, 4), t.getCellAsString(i, 27));
        }

        t = new RowTable (taxaGroup2FileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaGroupIDMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 7));
        }

        t = new RowTable (addInfileS);
        List<String> headerl = t.getHeader();

        try {
            BufferedWriter bw = IOUtils.getTextWriter(addOutfileS);
            for (int i = 0; i < headerl.size(); i++) {
                bw.write(headerl.get(i) + "\t");
            }
            bw.write("Group\tSubspecies\tGroupID\tRatio");
//            bw.write("Taxa\tDeleteriousCountPerHaplotype\tSiteCountWithMinDepth\tGroup\tSubspecies\tGroupID\tRatio");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) { //这里的t是 addInfileS
                int index = Arrays.binarySearch(taxa, t.getCellAsString(i, 0));
                if (index < 0) continue;
                if(taxaGroupMap.get(taxa[index]).equals("ExclusionHexaploid") || taxaGroupMap.get(taxa[index]).equals("ExclusionTetraploid")) continue;
                double genotypesite = Double.valueOf(t.getCellAsDouble(i, 2));
                if(genotypesite == 0) continue;
                double ratio = Double.valueOf(t.getCellAsDouble(i, 2))/Double.valueOf(t.getCellAsDouble(i, 3));
                for (int j = 0; j < t.getColumnNumber(); j++) {
                    bw.write(t.getCellAsString(i,j) + "\t");
                }
                bw.write(taxaGroupMap.get(taxa[index])
                        +"\t"+taxaSubMap.get(taxa[index])
                        +"\t"+taxaGroupIDMap.get(taxa[index])
                        +"\t"+String.format("%.3f",ratio));
                bw.newLine();
            }
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * 新的数据库分析，添加了dirivedSIFT的列信息,并且将每条染色体的情况都统计下来
     */
    public void countDeleteriousVMapII_byChr() {
        String delVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/002_exonSNPVCF"; //有害变异的VCF文件路径


//        String delVCFDirS = ""; //有害变异的VCF文件路径

//        String deleFileS = ""; //有害变异信息库
//        String addCountFileS = ""; //有害变异加性模型输出文件
//        String recCountFileS = ""; //有害变异隐形模型输出文件

//        String delVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/002_vcf/001_/del"; //有害变异的VCF文件路径
//        String deleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/001_/del_merged/chr_SNP_anno_del.txt.gz"; //有害变异信息库
//        String addCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/001_/additiveDeleterious_vmap2_bychr.txt"; //有害变异加性模型输出文件
//        String recCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/001_/recessiveDeleterious_vmap2_bychr.txt"; //有害变异隐形模型输出文件

//        String delVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/002_vcf/002_/del"; //有害变异的VCF文件路径
//        String deleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/002_/del_merged/chr_SNP_anno_del.txt.gz"; //有害变异信息库
//        String addCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/002_/additiveDeleterious_vmap2_bychr.txt"; //有害变异加性模型输出文件
//        String recCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/002_/recessiveDeleterious_vmap2_bychr.txt"; //有害变异隐形模型输出文件

//        String delVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/002_vcf/003_/syn"; //有害变异的VCF文件路径
//        String deleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/003_/syn_merged/chr_SNP_anno_syn.txt.gz"; //有害变异信息库
//        String addCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/003_/additiveDeleterious_vmap2_bychr.txt"; //有害变异加性模型输出文件
//        String recCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/003_/recessiveDeleterious_vmap2_bychr.txt"; //有害变异隐形模型输出文件

//        String delVCFDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/002_vcf/000_/del"; //有害变异的VCF文件路径
//        String deleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/001_snp/000_/del_merged/chr_SNP_anno_del.txt.gz"; //有害变异信息库
//        String addCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/000_/additiveDeleterious_vmap2_bychr.txt"; //有害变异加性模型输出文件
//        String recCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/024_deleteriousBiology/003_VMap2.1DelCount/001_delCount/000_/recessiveDeleterious_vmap2_bychr.txt"; //有害变异隐形模型输出文件




//        String deleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/000_/del_merged/chr_SNP_anno_del.txt.gz"; //有害变异信息库
//        String addCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/001_delCount/000_/additiveDeleterious_vmap2_bychr.txt"; //有害变异加性模型输出文件

//        String deleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/001_/del_merged/chr_SNP_anno_del.txt.gz"; //有害变异信息库
//        String addCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/001_delCount/001_/additiveDeleterious_vmap2_bychr.txt"; //有害变异加性模型输出文件

//        String deleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/002_/del_merged/chr_SNP_anno_del.txt.gz"; //有害变异信息库
//        String addCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/001_delCount/002_/additiveDeleterious_vmap2_bychr.txt"; //有害变异加性模型输出文件

//        String deleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/003_/syn_merge/chr_SNP_anno_syn.txt.gz"; //有害变异信息库
//        String addCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/001_delCount/003_/additiveDeleterious_vmap2_bychr.txt"; //有害变异加性模型输出文件

        String deleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/004_gerp/del_merged/chr_SNP_anno_del.txt.gz"; //有害变异信息库
        String addCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/001_delCount/004_gerp/additiveDeleterious_vmap2_bychr.txt"; //有害变异加性模型输出文件

//        String deleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/005_nonsynonymousGerp1/del_merged/chr_SNP_anno_del.txt.gz"; //有害变异信息库
//        String addCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/001_delCount/005_nonsynonymousGerp1/additiveDeleterious_vmap2_bychr.txt"; //有害变异加性模型输出文件

//        String deleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/006_synonymousGerp1/syn_merge/chr_SNP_anno_syn.txt.gz"; //有害变异信息库
//        String addCountFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/005_VMap2.1DelCount/001_delCount/006_synonymousGerp1/additiveDeleterious_vmap2_bychr.txt"; //有害变异加性模型输出文件


        int minDepth = 2;//inclusive
        int chrNum = 42;
        ArrayList<Integer> chrList = new ArrayList();
        for (int i = 0; i < chrNum; i++) {
            chrList.add(i + 1);
        }


        RowTable<String> t = new RowTable(deleFileS);
        TIntArrayList[] posList = new TIntArrayList[chrNum];
        TCharArrayList[] charList = new TCharArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) { //集合类数组，要初始化每一个list
            posList[i] = new TIntArrayList();
            charList[i] = new TCharArrayList();
        }
        new AoFile().readheader(deleFileS);
        String derivedAllele = null;
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = t.getCellAsInteger(i, 1) - 1; //染色体号的索引

            //################### 需要修改 //###################//###################//###################//###################
//            String ancestralAllele = t.getCell(i, 22); //不同的数据库，这一列的信息不一样，千万要注意!!!!!!!!!!!!!!!!! 祖先基因的数据库
            String ancestralAllele = t.getCell(i, 23); //不同的数据库，这一列的信息不一样，千万要注意!!!!!!!!!!!!!!!!! 祖先基因的数据库
            //################### 需要修改 //###################//###################//###################//###################

            String majorAllele = t.getCell(i, 5);
            String minorAllele = t.getCell(i, 6);
            String pos = t.getCell(i, 2);
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
            if (!ancestralAllele.equals(minorAllele) && !ancestralAllele.equals(majorAllele) && !ancestralAllele.equals("NA") && !ancestralAllele.equals("null")) {
                System.out.println(ancestralAllele+ "\t" + majorAllele + "\t" + minorAllele+ "\t" + index + "\t" + pos+ "\t" + (i+1));
            }
        }
        int[][] delePos = new int[chrNum][];
        char[][] deleChar = new char[chrNum][];
        for (int i = 0; i < chrNum; i++) { //将每一个list转化为数组
            delePos[i] = posList[i].toArray();
            deleChar[i] = charList[i].toArray();
            Arrays.sort(delePos[i]);
        }


        String vmap2TaxaList = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/001_taxaList/004_vmap2.1_TaxaList/TaxaList.txt";
        List<String> taxaList = new AoFile().getStringListwithoutHeader(vmap2TaxaList, 0);
        String[] taxa = taxaList.toArray(new String[taxaList.size()]);
        Arrays.sort(taxa);


        int taxaNum = taxa.length;
        double[][] addCount = new double[chrNum][taxa.length];
        int[][] recCount = new int[chrNum][taxa.length];
        int[][] siteWithMinDepthCount = new int[chrNum][taxa.length]; //每个taxa的深度

        chrList.parallelStream().forEach(chr -> {
            String delVmapFileS = "chr" + PStringUtils.getNDigitNumber(3, chr) + "_exon_vmap2.1.vcf.gz";
//            String delVmapFileS = "chr" + PStringUtils.getNDigitNumber(3, chr) + "_syn_vmap2.1.vcf.gz"; /////////此处需要修改

            delVmapFileS = new File(delVCFDirS, delVmapFileS).getAbsolutePath();
            BufferedReader br = IOUtils.getTextGzipReader(delVmapFileS);
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
                        int pos = Integer.valueOf(l.get(1));
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
                            //如果ref是derived allele，那么idx[0]=0;idx[1]=1. 当时0/0时，sum=2; 0/1时，sum=1; 1/1时， sum=0.
                            //如果alt是derived allele，那么idx[0]=1;idx[1]=0. 当时0/0时，sum=0; 0/1时，sum=1; 1/1时， sum=2.
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
            bw.write("Taxa\tChr\tDeleteriousCountPerHaplotype\tSiteCountWithMinDepth"); //每个taxa有多少个加性效应的derivedAllele 每个taxa在del库中含有基因型的个数
            bw.newLine();
            for (int i = 0; i < addCount.length; i++) {
                String chr = String.valueOf(i + 1);
                for (int j = 0; j < addCount[0].length; j++) {
                    bw.write(taxa[j] + "\t" + chr + "\t" + String.valueOf(addCount[i][j]) + "\t" + String.valueOf(siteWithMinDepthCount[i][j]));
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void mergeDelgenicSNPAnnotation() {

        //change
//        String infileDirS = "";
//        String outfileS = "";
//        new CountSites().mergeTxtbysuffix(infileDirS, outfileS, ".txt.gz");



//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/000_/del";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/000_/del_merged/chr_SNP_anno_del.txt.gz";
//        new CountSites().mergeTxtbysuffix(infileDirS, outfileS, ".txt.gz");

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/001_/del";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/001_/del_merged/chr_SNP_anno_del.txt.gz";
//        new CountSites().mergeTxtbysuffix(infileDirS, outfileS, ".txt.gz");

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/002_/del";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/002_/del_merged/chr_SNP_anno_del.txt.gz";
//        new CountSites().mergeTxtbysuffix(infileDirS, outfileS, ".txt.gz");

//                String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/003_/syn";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/003_/syn_merge/chr_SNP_anno_syn.txt.gz";
//        new CountSites().mergeTxtbysuffix(infileDirS, outfileS, ".txt.gz");

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/004_gerp/del";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/004_gerp/del_merged/chr_SNP_anno_del.txt.gz";
        new CountSites().mergeTxtbysuffix(infileDirS, outfileS, ".txt.gz");


    }


    public void addAnc(){

        //no change
//        String anceDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/003_ancestral/hv_brdis/002_byChrID";
//        String anceDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/003_ancestral/hv_br_jap/002_byChrID";
//        String anceDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/006_ancestralfromLipeng/002_byChrID";
        String anceDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/007_ancestral_Barley_secale/002_byChrID";

                //change
//        String infileDirS = "";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_deleteriousISite/001_syn";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_deleteriousISite/002_nonsyn";
//        String infileDirS ="/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_deleteriousISite/004_nonsynSiftDerived";
//                String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_deleteriousISite/003_nonsynSiftRefAlt";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_deleteriousISite/005_gerp";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/000_/del";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/001_/del";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/002_/del";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/003_/syn";
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/004_snp/004_gerp/del";

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/005_exonSNPAnnotation";
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/007_exonSNPAnnotation_addDAF_barleyUratu";
        List<File> fs =IOUtils.getVisibleFileListInDir(infileDirS);
        for (int i = 0; i < fs.size(); i++) {
            String infileS = fs.get(i).getAbsolutePath();
            String chr = new File(infileS).getName().substring(3,6);
//            String anceS = new File(anceDirS,"chr"+chr+".hovul.brdis.exon.probs.ancestral.txt.gz").getAbsolutePath();
//            String anceS = new File(anceDirS,"chr" + chr + ".hovul.brdis.orjap.exon.probs.ancestral.txt.gz").getAbsolutePath();
//            String anceS = new File(anceDirS,"chr" + chr + "_barleyVSsecale_ancestralAllele.txt.gz").getAbsolutePath();
            String anceS = new File(anceDirS,"chr" + chr + "_secer_hv_ancestral.txt.gz").getAbsolutePath();
            HashMap<Integer,String> hm = AoFile.getHashMapintKey(anceS,1,2);
            new AoFile().addColumbyint(infileS,2,hm,"Ancestral_barley_secale_parsimony");
        }
    }


    public void sampleOutdata(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/001_est-sfs";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/107_estsfs/002_sample";
        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        double extractRatio = 0.01;
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".probs")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().split(".probs")[0] + "_subset.txt.gz").getAbsolutePath();
                } else if (infileS.endsWith(".probs.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().split(".probs.gz")[0] + "_subset.txt.gz").getAbsolutePath();
                }
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                bw.write("MajorAlleleP\tMajorP");
                bw.newLine();
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    List<String> l = PStringUtils.fastSplit(temp);
                    if (temp.startsWith("0")) continue;
                    double r = Math.random();
                    double ratio = extractRatio;
                    if (r > ratio) {
                        continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                    }
                    String value = l.get(6);
                    List<String> ll = PStringUtils.fastSplit(value," ");
                    String v1 = ll.get(2);
                    double v2 = Double.parseDouble(ll.get(3));
                    double v3 = Double.parseDouble(ll.get(4));
                    double v4 = Double.parseDouble(ll.get(5));
                    double v5 = Double.parseDouble(ll.get(6));
                    double[] v = {v2,v3,v4,v5};
                    double max = Arrays.stream(v).max().getAsDouble();

                    bw.write(v1 + "\t" + String.format("%.3f",max));
                    bw.newLine();
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

}
