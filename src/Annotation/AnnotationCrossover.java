package Annotation;


import analysis.wheat.VMap2.VMapDBUtils;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;


import java.io.BufferedWriter;
import java.io.File;
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

    }

    /**
     *
     *
     */
    public void addRecombination () {
        String dirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/test/001/001_exonSNPAnnotation";
        String outDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/test/001/002_addRecombination";
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
        List<File> fList = IOUtils.getFileListInDirEndsWith(dirS, ".txt");
        fList.parallelStream().forEach(f -> {
            //String outfileS = new File (outDirS, f.getName()).getAbsolutePath();
            Dyad<String, List<String>> two = VMapDBUtils.getDBInfo(f.getAbsolutePath()); //返回Dyad类型
            String header = two.getFirstElement(); //返回表头
            List<String> recordList = two.getSecondElement(); //返回每一行的内容的集合
            String[] tem = header.split("\t");
            try {
                BufferedWriter bw = IOUtils.getTextGzipWriter(f.getAbsolutePath());
                StringBuilder sb = new StringBuilder(header);
                sb.append("\t").append("RecombinationRate");
                bw.write(sb.toString());
                bw.newLine();
                int chrIndex = -1;
                int posIndex = -1;
                int currentPos = -1;
                List<String> l  = null;
                for (int i = 0; i < recordList.size(); i++) {
                    sb.setLength(0);
                    l = PStringUtils.fastSplit(recordList.get(i)); //读每一行的内容
                    chrIndex = Integer.parseInt(l.get(1))-1; // 索引1 含有染色体号
                    currentPos = Integer.parseInt(l.get(2)); //索引2 含有位置
                    posIndex = startLists[chrIndex].binarySearch(currentPos); //在 刚刚的起始集合里搜索 index
                    if (posIndex < 0) posIndex = -posIndex-2;
                    if (posIndex < 0) {  //如果index小于0，说明在最开始区间的前面（第一个window的最前面），即没有在起始位点集合中
                        sb.append(recordList.get(i)).append("\t").append("NA");
                        bw.write(sb.toString());
                        bw.newLine();
                        continue;
                    }
                    if (currentPos < endLists[chrIndex].get(posIndex)) { //改点和所在window的右侧最大值相比，如果小于他，说明是在该window里
                        sb.append(recordList.get(i)).append("\t").append(crossLists[chrIndex].get(posIndex)); //获取该窗口对应的 重组率值
                    }
                    else {
                        sb.append(recordList.get(i)).append("\t").append("NA"); //如果该点不在window的右侧范围，说明是位于最后一个window并且超出最后一个window的右侧最大值
                    }
                    bw.write(sb.toString());
                    bw.newLine();
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
