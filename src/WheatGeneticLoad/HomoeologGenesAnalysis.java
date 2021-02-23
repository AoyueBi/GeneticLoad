package WheatGeneticLoad;

import AoUtils.AoColor;
import AoUtils.AoFile;
import AoUtils.AoMath;
import gnu.trove.list.array.TDoubleArrayList;
import pgl.infra.utils.PStringUtils;

import java.awt.image.Kernel;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class HomoeologGenesAnalysis {

    public HomoeologGenesAnalysis(){
        this.getAverageDistance();
//        this.getCol();
//        this.extractDistance();




    }

    public void extractDistance(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/007_Global/002_hexaploid_byIndividual/002_addAvedistance";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/007_Global/002_hexaploid_byIndividual/003_extractAve";
        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        AoFile.readheader(fsList.get(0).getAbsolutePath());
        int a=3;
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS, f.getName()).getAbsolutePath();
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                bw.write("AveDistance");
                bw.newLine();
                String header = br.readLine();
                String temp = null;
                List<String> l = new ArrayList<>();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    for (int i = 0; i < 5; i++) {


                    }
                    String ave = l.get(2010);
                    bw.write(ave);
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


    public void getCol(){
        String [] key = {"M011","M110"};
        AoColor.sort(key);
        AoColor.models7(key);

    }

    /**
     * 获取每个群体到gloable的平均距离
     * 先确定那几列的值需要进行计算，再建立 DoubleList,最后求平均值
     *
     */
    public void getAverageDistance(){
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/007_Global/001_hexaploid_perCDSperGenotype/001_source";
//        String outfileDirs = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/007_Global/001_hexaploid_perCDSperGenotype/002_addAvedistance";

        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/007_Global/002_hexaploid_byIndividual/001_source";
        String outfileDirs = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/108_geneDB/007_Global/002_hexaploid_byIndividual/002_addAvedistance";
        new File(outfileDirs).mkdirs();

        List<File> fsList = AoFile.getFileListInDir(infileDirS);

//        String testFileS = fsList.get(0).getAbsolutePath();
//        int cntColumn = AoFile.countFileColumnNumber(testFileS);
//        int cntList = cntColumn/5-1; //有多少个群体数
//        int[] indexCal = new int[cntList];
//        for (int i = 0; i < cntList; i++) { //从第一个群体开始计算index，依次类推
//            indexCal[i] = (i+2)*5-1;
//        }
//        Arrays.sort(indexCal);
//        int a=3;
        fsList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = new File(outfileDirs,f.getName()).getAbsolutePath();
            try {
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String temp = null;
                String header = br.readLine();
                bw.write(header + "\tAveDistance");
                bw.newLine();
                List<String> l = new ArrayList<>();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    cnt++;
                    TDoubleArrayList disList = new TDoubleArrayList();
                    for (int i = 0; i < l.size(); i++) { //注意，在达兴的程序里，任何群体为NA和M000模型的都不计算距离
                        if (i<5)continue;
                        if ((i+1)%5 == 0){
                            String value = l.get(i);
                            if (value.startsWith("N"))continue;
                            disList.add(Double.parseDouble(l.get(i)));
                        }
                    }
                    String ave = AoMath.getRelativeMean(disList);
                    bw.write(temp + "\t" + ave);
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
        });
    }
}
