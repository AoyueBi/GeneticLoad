package WheatGeneticLoad;

import AoUtils.AoFile;
import gnu.trove.list.array.TDoubleArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import xiaohan.analysis.RNAseq.NewFile;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class RebuildVCF {

    public RebuildVCF(){
//        this.statVcfDepth_SD_PValue();
        this.mergeTxt();
    }

    public void mergeTxt(){
        AoFile.mergeTxt("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/002_depthCal/","/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/001_merge/001_vmap2_subset0.001_depth.txt.gz");

    }

/**
 * 这里采用单线程 对抽样*的vcf文件进行每个位点深度和sd进行统计和Pvalue获取，制成一个表格，chr pos averageDepth SD
 * PValue 注意表格不要以#开头，否则被注释，看不到

 */

    public void statVcfDepth_SD_PValue() {

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/001_subsetData";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/002_depthCal";
        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f  -> {
            String infileS = f.getAbsolutePath();
            String outfileS = new File(outfileDirS,f.getName().split(".vcf")[0] + "_depth.txt.gz").getAbsolutePath();

            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String[] taxa = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            try {
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("##"))continue;
                    if(temp.startsWith("#C")){
                        List<String> linetaxa = PStringUtils.fastSplit(temp, "\t");
                        bw.write(linetaxa.get(0).replaceFirst("#", "") + "\t" + linetaxa.get(1) + "\t" + "AverageDepth\tSD");
                        bw.newLine();
                        taxa = new String[linetaxa.size() - 9];
                        continue;
                    }
                    cnt++; // 对snp开始计数
                    if (cnt % 1000000 == 0) {
                        System.out.println(String.valueOf(cnt) + " lines");
                    }
                    l = PStringUtils.fastSplit(temp, "\t");
                    String chr = l.get(0);
                    String pos = l.get(1);
                    TDoubleArrayList depthList = new TDoubleArrayList();
                    for (int i = 0; i < taxa.length; i++) {
                        String genoS = l.get(i + 9);
                        if (genoS.startsWith(".")) {
                            depthList.add(0);
                            continue;
                        }
                        List<String> ll = PStringUtils.fastSplit(genoS, ":");
                        List<String> lll = PStringUtils.fastSplit(ll.get(1), ",");
                        int depth = Integer.valueOf(lll.get(0)) + Integer.valueOf(lll.get(1));
                        depthList.add(depth);
                    }
                    double[] dep = depthList.toArray();
                    DescriptiveStatistics d = new DescriptiveStatistics(dep);
                    double relativeMean = d.getMean();
                    double sd = d.getStandardDeviation();
                    //计算完毕，接下来开始写入文件
                    StringBuilder sb = new StringBuilder();
                    sb.append(chr).append("\t").append(pos).append("\t").append(String.format("%.6f", relativeMean)).append("\t").append(String.format("%.6f", sd));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                System.out.println(temp);
                e.printStackTrace();
                System.exit(1);
            }
            System.out.println(infileS + " is calculated well done");

        });

    }
}
