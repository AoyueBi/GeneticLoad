package PopulationAnalysis;

import AoUtils.AoFile;
import analysis.wheatVMap2.VMapDBUtils;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import pgl.format.table.ColumnTable;
import pgl.format.table.RowTable;
import pgl.utils.Dyad;
import pgl.utils.IOUtils;
import pgl.utils.PStringUtils;
import pgl.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class XPCLR {
    public XPCLR(){
//        this.convertCoordinate();
//        this.test();
//        this.addgeneticPos();
//        this.calDensity("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/002_snp/chr001.subgenome.maf0.01.SNP_bi.subset.pos.Base.txt.gz",1,3,2000000,2000000,"/Users/Aoyue/Documents/test.txt");

        this.getGenotypeXPCLR();
    }

    /**
     *
     *
     */
    public void getGenotypeXPCLR(){
        String popfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/003_genoTest/001_pop/Cultivar.txt";
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/022_subsetVCF/001_singleChr0.001/chr001.subgenome.maf0.01.SNP_bi.subset.vcf.gz";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/003_genoTest/002_geno/chr001_cultivar.geno.txt";

        List<String> queryTaxal = new AoFile().getStringListwithoutHeader(popfileS,0); //获取pop列表
        List<Integer> indexHexa = new ArrayList<>();

        try {
            BufferedReader br = null;
            if (infileS.endsWith(".vcf")) {
                br = IOUtils.getTextReader(infileS);
            } else if (infileS.endsWith(".vcf.gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            BufferedWriter bw = null;
            if (outfileS.endsWith(".geno.txt")){
                bw=IOUtils.getTextWriter(outfileS);
            }else if (outfileS.endsWith(".geno.txt.gz")) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }

            String temp = null;
            List<String> l = new ArrayList<>();
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
            BufferedWriter bw = null;
            if (outfileS.endsWith(".txt")){
                bw = IOUtils.getTextWriter(outfileS);
            }
            if (outfileS.endsWith(".txt.gz")){
                bw = IOUtils.getTextGzipWriter(outfileS);
            }
            bw.write("CHROM\tBIN_START\tBIN_END\tN_VARIANTS");
            bw.newLine();
            for (int i = 0; i < bound.length; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(chr).append("\t").append(bound[i][0]).append("\t").append(bound[i][1]-1).append("\t").append(count[i]);
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
        while (start < chrLength) {
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
     * Today I  just read REF and not programming
     * Today I  just made annual summary and not programming
     */

    public void test(){
        int a =12; int b = 12;
        int c = Math.max(a,b);
        float d = 0.131f;
        System.out.println(String.format("%.6f", d/100));
        System.out.println(d/100);

    }

    /**
     *
     * 根据science 发表的重组文件，向 snp info 中添加 genetic position
     */
    public void addgeneticPos () {
//        String dirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/019_popGen/104_XPCLR/002_snp";
//        String recombinationFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/007_recombination/001_recombination/iwgsc_refseqv1.0_mapping_data_chrID.txt";

        String dirS = "/data4/home/aoyue/vmap2/analysis/022_XPCLR/001_prepare/001_chrposRefAlt";
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

}
