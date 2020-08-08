package Annotation;


import AoUtils.AoFile;
import AoUtils.WheatUtils;
import analysis.wheat.VMap2.VMapDBUtils;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import pgl.infra.table.ColumnTable;
import pgl.infra.table.RowTable;
import pgl.infra.utils.Dyad;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;


/**
 * @author AoyueBi
 * @data 2020-05-31 23:01
 */
public class AnnotationCrossover {
    public AnnotationCrossover(){
//        this.convertCoordinate();
        this.addRecombination();
//        this.mergeExonAnnotation();

    }

    /**
     * 获取重组率在标准化的坐标上的分布
     */
    public void getScaledPos(){
        String infileS = "/Users/Aoyue/Documents/Data/wheat/article/iwgsc_refseqv1.0_recombination_rate_analysis/iwgsc_refseqv1.0_recombination_rate.txt";
        String outfileS = "/Users/Aoyue/Documents/Data/wheat/article/iwgsc_refseqv1.0_recombination_rate_analysis/iwgsc_refseqv1.0_recombination_rate_addScalePos.txt";
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String header = br.readLine();
            bw.write(header + "\tPosScale");
            bw.newLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;

                String chromosome = l.get(0).substring(3);
                int posOnchromosome = Integer.parseInt(l.get(1));
                String posScaled = WheatUtils.getScaledPos(chromosome,posOnchromosome);
                bw.write(temp + "\t" + posScaled);
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

    public void mergeExonAnnotation(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/020_exonSNPAnnotation_addRecombination_from018";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/021_exonSNPAnnotation_addRecombination_from018_merge/001_exonSNP_anno.txt.gz";
        AoFile.mergeTxt(infileDirS,outfileS);
    }

    /**
     *
     * 在 exonAnnotation 数据库中添加 recombination 一列
     */
    public void addRecombination () {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/018_exonSNPAnnotation";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/genicSNP/020_exonSNPAnnotation_addRecombination_from018";
        String recombinationFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/007_recombination/001_recombination/iwgsc_refseqv1.0_recombination_rate_chrID.txt";
        ColumnTable<String> t = new ColumnTable<>(recombinationFileS);
        int chrNum = Integer.parseInt(t.getCell(t.getRowNumber()-1, 0)); //获取最后一行第0列的数字，即染色体最大值，这里是42号染色体
        TIntArrayList[] startLists = new TIntArrayList[chrNum]; //42条染色体中，每条染色体的窗口起始位置的集合
        TIntArrayList[] endLists = new TIntArrayList[chrNum]; //42条染色体中，每条染色体的窗口结束位置的集合
        TFloatArrayList[] crossLists = new TFloatArrayList[chrNum]; //42条染色体中，每条染色体的每个窗口对应cross数值集合
        for (int i = 0; i < startLists.length; i++) { //对每个数组内的集合进行初始化
            startLists[i] = new TIntArrayList();
            endLists[i] = new TIntArrayList();
            crossLists[i] = new TFloatArrayList();
        }
        int index = -1;
        for (int i = 0; i < t.getRowNumber(); i++) {
            index = Integer.parseInt(t.getCell(i, 0))-1; //染色体号的索引，即1号染色体索引为0
            startLists[index].add(Integer.parseInt(t.getCell(i, 1))); //将每个Bin的起始位置加入集合中
            endLists[index].add(Integer.parseInt(t.getCell(i, 2)));
            crossLists[index].add(Float.parseFloat(t.getCell(i, 3)));
        }
        List<File> fList = AoFile.getFileListInDir(infileDirS);
        fList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = new File(outfileDirS,f.getName()).getAbsolutePath();
            try {
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String header = br.readLine();
                bw.write(header + "\tRecombinationRate");bw.newLine();
                String temp = null;
                List<String> l = new ArrayList<>();
                int cnt = 0;
                int chrIndex = -1;
                int posIndex = -1;
                int currentPos = -1;
                StringBuilder sb = new StringBuilder();
                while ((temp = br.readLine()) != null) {
                    sb.setLength(0);
                    l = PStringUtils.fastSplit(temp);
                    chrIndex = Integer.parseInt(l.get(1))-1; // 索引1 含有染色体号
                    currentPos = Integer.parseInt(l.get(2)); //索引2 含有位置
                    posIndex = startLists[chrIndex].binarySearch(currentPos); //在 刚刚的起始集合里搜索 index
                    if (posIndex < 0) posIndex = -posIndex-2;
                    if (posIndex < 0) {  //如果index小于0，说明在最开始区间的前面（第一个window的最前面），即没有在起始位点集合中
                        sb.append(temp).append("\t").append("NA");
                        bw.write(sb.toString());
                        bw.newLine();
                        continue;
                    }
                    if (currentPos < endLists[chrIndex].get(posIndex)) { //该点和所在window的右侧最大值相比，如果小于他，说明是在该window里
                        sb.append(temp).append("\t").append(crossLists[chrIndex].get(posIndex)); //获取该窗口对应的 重组率值
                    }
                    else {
                        sb.append(temp).append("\t").append("NA"); //如果该点不在window的右侧范围，说明是位于最后一个window并且超出最后一个window的右侧最大值
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                    cnt++;
                }
                br.close();
                bw.flush();
                bw.close();
                System.out.println(f.getName() + "\tis completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    /**
     * 通过 Chr1A- Chr7D pos 信息转换为 chrID 格式
     *
     */
    public void convertCoordinate () {
        String infileS = "/Users/Aoyue/Documents/Data/wheat/article/iwgsc_refseqv1.0_recombination_rate_analysis/iwgsc_refseqv1.0_recombination_rate.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/007_recombination/001_recombination/iwgsc_refseqv1.0_crossover_rate_chrID.txt";
        RowTable<String> t = new RowTable<>(infileS);
        String header = "ChrID\tPosStart\tPosEnd\tCrossover_frequency";
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(header);
            bw.newLine();
            String chromosome = null;
            int start = -1;
            int end = -1;
            int chr1 = -1;
            int chr2 = -1;
            int pos1 = -1;
            int pos2 = -1;
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < t.getRowNumber(); i++) {
                sb.setLength(0);
                chromosome = t.getCell(i,0).replaceFirst("chr", "");
                start = Integer.parseInt(t.getCell(i,1));
                end = Integer.parseInt(t.getCell(i,2));
                chr1 = RefV1Utils.getChrID(chromosome, start);
                chr2 = RefV1Utils.getChrID(chromosome, end);
                if (chr1 == chr2) { //如果所在区间的起始和末尾都属于一条染色体
                    pos1 = RefV1Utils.getPosOnChrID(chromosome, start);
                    pos2 = RefV1Utils.getPosOnChrID(chromosome, end);
                    sb.append(chr1).append("\t").append(pos1).append("\t").append(pos2).append("\t").append(t.getCell(i, 4));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                else {
                    pos1 = RefV1Utils.getPosOnChrID(chromosome, start);
                    sb.append(chr1).append("\t").append(pos1).append("\t").append(RefV1Utils.getChrIDLength(chr1)+1).append("\t").append(t.getCell(i, 4));
                    bw.write(sb.toString());
                    bw.newLine();
                    sb.setLength(0);
                    pos2 = RefV1Utils.getPosOnChrID(chromosome, end);
                    sb.append(chr2).append("\t").append(1).append("\t").append(pos2).append("\t").append(t.getCell(i, 4));
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
