package WheatGeneticLoad;

import AoUtils.AoFile;
import gnu.trove.list.array.TDoubleArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.graph.r.DensityPlot;
import pgl.graph.r.ScatterPlot;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

public class RebuildVCF {

    public RebuildVCF(){
//        this.statVcfDepth_SD_PValue();
//        this.mergeTxt();
        /**
         * Fei's test
         */
//        this.densityFilter();
//        this.mkReliableGenotypeSite();
//        this.mkReliableIntersection();
//        this.intersectionCheck();


//        this.addReliableGroupfromSample();
//        this.addReliableInsecterGroupfromSample();
//        this.mergeTxt2();

        /**
         * 对fastcall 文件的检查
         */

//        this.checkErrorFromFastCall();

    }

    /**
     * 在对fastcall文件进行VCF提取时，发现index out of bound, 故写程序提取结果查找异常。
     */
    public void checkErrorFromFastCall(String infileS){
        //model
//        String infileS = "";

//        String infileS = "/Users/Aoyue/Documents/test.vcf";
        try {
            BufferedReader br = AoFile.readFile(infileS);
            String header = br.readLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            int cntError = 0;
            int size = Integer.MIN_VALUE;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##"))continue;
                if (temp.startsWith("#C")){
                    l = PStringUtils.fastSplit(temp);
                    size = l.size();
                    System.out.println("total size is " + size);
                    continue;
                }
                cnt++;
                l = PStringUtils.fastSplit(temp);
                if (l.size() == size)continue;
                else {
                    System.out.println(temp);
                    cntError++;
                }
            }
            br.close();
            System.out.println(cnt + "  snps    " + cntError + "    error snps in " + new File(infileS).getName());
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        //java -jar PlantGenetics.jar /data4/home/aoyue/vmap2/genotype/fastcall/abd/chr002.ABDgenome.vcf.gz > log_20201111/log_chr002_10X_fastcall_hexaploid.txt 2>&1 &
        //java -jar PlantGenetics.jar /data4/home/aoyue/vmap2/genotype/fastcall/abd/chr033.ABDgenome.vcf.gz > log_20201111/log_chr033_10X_fastcall_hexaploid.txt 2>&1 &
        //java -jar PlantGenetics.jar /data4/home/aoyue/vmap2/genotype/fastcall/abd/chr034.ABDgenome.vcf.gz > log_20201111/log_chr034_10X_fastcall_hexaploid.txt 2>&1 &
        // sed -i '3353399d' chr002.ABDgenome.vcf &

        // 2 3353399
        //33 12060211
        //34 9975515

        //
    }





    public void mergeTxt2(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/004_AddReliableIntersectGroup";
        File[] fs = new File(infileDirS).listFiles();

//        fs = IOUtils.listFilesEndsWith(fs,"_AB_sample.txt.gz");
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/005_changeChrPos/AB_Popdepth_sample.txt.gz";

//        fs = IOUtils.listFilesEndsWith(fs,"_ABD_sample.txt.gz");
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/005_changeChrPos/ABD_Popdepth_sample.txt.gz";

        fs = IOUtils.listFilesEndsWith(fs,"_D_sample.txt.gz");
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/005_changeChrPos/D_Popdepth_sample.txt.gz";

        try{
            Arrays.sort(fs);
            String infileS = fs[0].getAbsolutePath();
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            //read header
            bw.write("Chr\t" + br.readLine());
            bw.newLine();

            int cnttotal=0;
            for (int i = 0; i < fs.length; i++) {
                infileS = fs[i].getAbsolutePath();
                int chrID = Integer.parseInt(fs[i].getName().substring(3,6));
                br = AoFile.readFile(infileS);
                br.readLine();
                String temp = null; //read header
                int cnt = 0;
                List<String> l = new ArrayList<>();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    cnt++;
                    cnttotal++;
                    int pos = Integer.parseInt(l.get(0));
                    String chr = RefV1Utils.getChromosome(chrID,pos);
                    int posOnChrosome = RefV1Utils.getPosOnChromosome(chrID,pos);
                    bw.write(chr + "\t" + posOnChrosome);
                    for (int j = 1; j < l.size(); j++) {
                        bw.write("\t" + l.get(j));
                    }
                    bw.newLine();
                }
                System.out.println(fs[i].getName() + "\t" + cnt);
            }
            System.out.println("Total lines without header count is " + cnttotal + " at merged file " + outfileS );
            br.close();
            bw.flush();
            bw.close();

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     * 为抽样的文件加上分组信息，1 代表即属于 reliable 又属于 insertion 的位点， 0 代表不属于 intersection 的位点
     */
    public void addReliableInsecterGroupfromSample(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/003_AddReliableGroup";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/004_AddReliableIntersectGroup";
        String reliableDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/reliableSites/ABD_intersect";

        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String chr = f.getName().split("_")[0];
            String reliableFileS = new File(reliableDirS,chr + "_intersect_reliable.txt.gz").getAbsolutePath();
            String outfileS = new File(outfileDirS,f.getName()).getAbsolutePath();
            this.intersection(infileS,reliableFileS,outfileS);
        });
    }

    private void intersection (String infileS, String reliableFile, String outfileS) {
        RowTable<String> t = new RowTable<>(infileS);
        int[] positions = t.getColumnAsIntArray(0);
        boolean[] ifSelected = new boolean[positions.length];
        try {
            BufferedReader br = IOUtils.getTextGzipReader(reliableFile);
            String temp = br.readLine();
            int currentPos = 0;
            int value;
            int index;

            while ((temp = br.readLine()) != null) {
                currentPos++;
                if (currentPos%10000000 == 0) System.out.println(currentPos);
                value = Integer.parseInt(temp);
                if (value == 0) {
                    continue;
                }
                index = Arrays.binarySearch(positions, currentPos);
                if (index < 0) {
                    continue;
                }
                ifSelected[index] = true;
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }

        try{
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String header = br.readLine();
            bw.write(header + "\tIfIntersect");
            bw.newLine();
            String temp = null;
            String ifIntersect = null;
            int cnt = -1;
            while ((temp = br.readLine()) != null) {
                cnt++;
                if (!ifSelected[cnt]){
                    ifIntersect="0";
                }
                if (ifSelected[cnt]){
                    ifIntersect="1";
                }
                bw.write(temp + "\t" + ifIntersect);
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();

        }
    }

    /**
     * 为抽样的文件加上分组信息，1 代表即属于 reliable 的位点， 0 代表不属于 reliable 的位点
     */

    public void addReliableGroupfromSample(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/002_split";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/003_AddReliableGroup";
        double proportionOfSite = 0.70;
        String sample1FileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/001_merged/abPopDep_sample.txt.gz";
        String sample2FileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/001_merged/abdPopDep_sample.txt.gz";
        String sample3FileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/001_merged/dPopDep_sample.txt.gz";

        File[] fs = AoFile.getFileArrayInDir(infileDirS);
        for (int i = 0; i < fs.length; i++) {
            String infileS = fs[i].getAbsolutePath();
            String outfileS = new File(outfileDirS,fs[i].getName()).getAbsolutePath();
            String ploid = fs[i].getName().split("_")[1];
//            int index = ploid.indexOf("Pop");
//            ploid = ploid.substring(0,index);
            if (ploid.equals("AB")){
                this.densityFilterAB(proportionOfSite,sample1FileS,infileS,outfileS);
            }else if (ploid.equals("ABD")){
                this.densityFilterABD(proportionOfSite,sample2FileS,infileS,outfileS);
            }else if (ploid.equals("D")){
                this.densityFilterD(proportionOfSite,sample3FileS,infileS,outfileS);
            }
        }
    }

    public void densityFilterD(double proportionOfSite,String sampleFileS, String infileS,String outfileS) {
        RowTable<String> t = new RowTable<>(sampleFileS);
        double depthStart = 3;
        double depthEnd = 17;
        double SDStart = 3;
        double SDEnd = 10;
        int nBin = 100;
        Grid gr = new Grid(depthStart, depthEnd, SDStart, SDEnd, nBin); //step1:new Grid
        for (int i = 0; i < t.getRowNumber(); i++) { //step2 : add site
            gr.addXY(t.getCellAsDouble(i,1), t.getCellAsDouble(i,2));
        }
        gr.buildHashMap(); //step3 : sort and build new hashmap based on cell count
        int indexThresh = gr.getOrderIndexOfProportionOfSite(proportionOfSite); //step4 : get index that in high density

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String header = br.readLine();
            bw.write(header + "\tIfReliable");
            bw.newLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                double x = Double.parseDouble(l.get(1));
                double y = Double.parseDouble(l.get(2));
                if (!gr.isHighDensity(x, y , indexThresh)){  //step5 : if high
                    bw.write(temp + "\t0");
                    bw.newLine();
                }else{
                    bw.write(temp + "\t1");
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

    public void densityFilterABD(double proportionOfSite,String sampleFileS, String infileS,String outfileS) {
        RowTable<String> t = new RowTable<>(sampleFileS);
        double depthStart = 3;
        double depthEnd = 15;
        double SDStart = 3;
        double SDEnd = 8;
        int nBin = 100;
        Grid gr = new Grid(depthStart, depthEnd, SDStart, SDEnd, nBin); //step1:new Grid
        for (int i = 0; i < t.getRowNumber(); i++) { //step2 : add site
            gr.addXY(t.getCellAsDouble(i,1), t.getCellAsDouble(i,2));
        }
        gr.buildHashMap(); //step3 : sort and build new hashmap based on cell count
        int indexThresh = gr.getOrderIndexOfProportionOfSite(proportionOfSite); //step4 : get index that in high density

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String header = br.readLine();
            bw.write(header + "\tIfReliable");
            bw.newLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                double x = Double.parseDouble(l.get(1));
                double y = Double.parseDouble(l.get(2));
                if (!gr.isHighDensity(x, y , indexThresh)){  //step5 : if high
                    bw.write(temp + "\t0");
                    bw.newLine();
                }else{
                    bw.write(temp + "\t1");
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

    /**
     * 自己写的
     * @param proportionOfSite
     * @param sampleFileS
     * @param infileS
     * @param outfileS
     */
    public void densityFilterAB(double proportionOfSite, String sampleFileS,String infileS,String outfileS) {
        RowTable<String> t = new RowTable<>(sampleFileS);
        double depthStart = 2;
        double depthEnd = 8;
        double SDStart = 2;
        double SDEnd = 8;
        int nBin = 100;
        Grid gr = new Grid(depthStart, depthEnd, SDStart, SDEnd, nBin); //step1:new Grid
        for (int i = 0; i < t.getRowNumber(); i++) { //step2 : add site
            gr.addXY(t.getCellAsDouble(i,1), t.getCellAsDouble(i,2));
        }
        gr.buildHashMap(); //step3 : sort and build new hashmap based on cell count
        int indexThresh = gr.getOrderIndexOfProportionOfSite(proportionOfSite); //step4 : get index that in high density

        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String header = br.readLine();
            bw.write(header + "\tIfReliable");
            bw.newLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;
                double x = Double.parseDouble(l.get(1));
                double y = Double.parseDouble(l.get(2));
                if (!gr.isHighDensity(x, y , indexThresh)){  //step5 : if high
                    bw.write(temp + "\t0");
                    bw.newLine();
                }else{
                    bw.write(temp + "\t1");
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


    /**
     * 以下是老师写的程序，以上都是自己写的程序
     */
    public void intersectionCheck () {
        String abSampleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/abPopDep_sample.txt";
        String abdSampleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/abdPopDep_sample.txt";
        String dSampleFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/dPopDep_sample.txt";

        String abFile1S = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/intersection/abPopDep_sample_chr001.txt";
        String abFile2S = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/intersection/abPopDep_sample_chr002.txt";
        String abdFile1S = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/intersection/abdPopDep_sample_chr001.txt";
        String abdFile2S = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/intersection/abdPopDep_sample_chr002.txt";
        String dFile1S = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/intersection/dPopDep_sample_chr005.txt";
        String dFile2S = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/intersection/dPopDep_sample_chr006.txt";

//        this.split(abSampleFileS, abFile1S, abFile2S);
//        this.split(abdSampleFileS, abdFile1S, abdFile2S);
//        this.split(dSampleFileS, dFile1S, dFile2S);

        String reliableChr001 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/001_reliableSites/ABD_intersect/chr001_intersect_reliable.txt.gz";
        String reliableChr002 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/001_reliableSites/ABD_intersect/chr002_intersect_reliable.txt.gz";
        String reliableChr005 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/001_reliableSites/ABD_intersect/chr005_intersect_reliable.txt.gz";
        String reliableChr006 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/001_depth/004_50000_Sites/001_reliableSites/ABD_intersect/chr006_intersect_reliable.txt.gz";

        String abDepthDensity = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/plot/intersection/abPopDep_intersect_depth_density.pdf";
        String abScatter = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/plot/intersection/abPopDep_intersect_scatter.pdf";
        String abdDepthDensity = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/plot/intersection/abdPopDep_intersect_depth_density.pdf";
        String abdScatter = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/plot/intersection/abdPopDep_intersect_scatter.pdf";
        String dDepthDensity = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/plot/intersection/dPopDep_intersect_depth_density.pdf";
        String dScatter = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/027_Rebuild_VMap2_VCF/plot/intersection/dPopDep_intersect_scatter.pdf";
//        this.plotIntersection(abFile1S, abFile2S, reliableChr001, reliableChr002, abDepthDensity, abScatter, "AB");
//        this.plotIntersection(abdFile1S, abdFile2S, reliableChr001, reliableChr002, abdDepthDensity, abdScatter, "ABD");
        this.plotIntersection(dFile1S, dFile2S, reliableChr005, reliableChr006, dDepthDensity, dScatter, "D");
    }

    private void plotIntersection (String file1, String file2, String reliableFile1, String reliableFile2, String densityPlot, String scatterPlot, String genomeType) {
        TDoubleArrayList depthList = new TDoubleArrayList();
        TDoubleArrayList sdList = new TDoubleArrayList();
        int totalCnt = 0;
        int depthFilterCnt = 0;
        int[] values =this.plotIntersection2(file1, reliableFile1, depthList, sdList);
        totalCnt+=values[0]; depthFilterCnt+=values[1];
        values = this.plotIntersection2(file2, reliableFile2, depthList, sdList);
        totalCnt+=values[0]; depthFilterCnt+=values[1];
        double[] x = depthList.toArray();
        double[] y = sdList.toArray();
        ScatterPlot s = new ScatterPlot(x, y);
        DensityPlot d = new DensityPlot(x);
        if (genomeType.equals("AB")) {
            s.setTitle("AB_intersect");
            s.setXLim(0, 8);
            s.setYLim(0, 8);
            s.setColor(255, 0, 0, 5);
            s.setXLab("Mean of depth");
            s.setYLab("SD of depth");
            s.saveGraph(scatterPlot);
            d.setXLim(0, 8);
            d.setXLab("Mean of depth");
            d.setYLab("Density");
            d.setTitle("AB_intersect");
            d.saveGraph(densityPlot);
        }
        else if (genomeType.equals("ABD")) {
            s.setTitle("ABD_intersect");
            s.setXLim(0, 20);
            s.setYLim(0, 12);
            s.setColor(255, 0, 0, 5);
            s.setXLab("Mean of depth");
            s.setYLab("SD of depth");
            s.saveGraph(scatterPlot);
            d.setXLim(0, 20);
            d.setXLab("Mean of depth");
            d.setYLab("Density");
            d.setTitle("ABD_intersect");
            d.saveGraph(densityPlot);
        }
        else if (genomeType.equals("D")) {
            s.setTitle("D_intersect");
            s.setXLim(0, 20);
            s.setYLim(0, 12);
            s.setColor(255, 0, 0, 5);
            s.setXLab("Mean of depth");
            s.setYLab("SD of depth");
            s.saveGraph(scatterPlot);
            d.setXLim(0, 20);
            d.setXLab("Mean of depth");
            d.setYLab("Density");
            d.setTitle("D_intersect");
            d.saveGraph(densityPlot);
        }
        System.out.println(genomeType);
        System.out.println((double)x.length/totalCnt);

    }

    private int[] plotIntersection2 (String file, String reliableFile, TDoubleArrayList depthList, TDoubleArrayList sdList) {
        RowTable<String> t = new RowTable<>(file);
        int[] positions = t.getColumnAsIntArray(0);
        boolean[] ifSelected = new boolean[positions.length];
        int depthFilterCnt = 0;
        try {
            BufferedReader br = IOUtils.getTextGzipReader(reliableFile);
            String temp = br.readLine();
            int currentPos = 0;
            int value;
            int index;

            while ((temp = br.readLine()) != null) {
                currentPos++;
                if (currentPos%10000000 == 0) System.out.println(currentPos);
                value = Integer.parseInt(temp);
                if (value == 0) {
                    continue;
                }
                index = Arrays.binarySearch(positions, currentPos);
                if (index < 0) {

                    continue;
                }
                ifSelected[index] = true;
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (!ifSelected[i]) continue;
            depthList.add(Double.parseDouble(t.getCell(i,1)));
            sdList.add(Double.parseDouble(t.getCell(i,2)));
        }
        int[] values = new int[2];
        values[0] = t.getRowNumber();
        values[1] = t.getRowNumber()-depthFilterCnt;
        return values;
    }

    private void split (String sourceFileS, String file1, String file2) {
        try {
            BufferedReader br = IOUtils.getTextReader(sourceFileS);
            String header  = br.readLine();
            String temp = null;
            BufferedWriter bw1 = IOUtils.getTextWriter(file1);
            BufferedWriter bw2 = IOUtils.getTextWriter(file2);
            bw1.write(header);
            bw1.newLine();
            bw2.write(header);
            bw2.newLine();
            List<String> l = new ArrayList<>();
            int current = -1;
            int next = -1;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                next = Integer.parseInt(l.get(0));
                if (next < current) {
                    bw2.write(temp);
                    bw2.newLine();
                    while ((temp = br.readLine()) != null) {
                        bw2.write(temp);
                        bw2.newLine();
                    }
                }
                else {
                    bw1.write(temp);
                    bw1.newLine();
                }
                current = next;
            }
            bw1.flush();
            bw1.close();
            bw2.flush();
            bw2.close();
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }


    public void mkReliableIntersection () {
        String abDirS = "/Volumes/VMap2_Fei/reliableSites/AB/";
        String abdDirS = "/Volumes/VMap2_Fei/reliableSites/ABD/";
        String dDirS = "/Volumes/VMap2_Fei/reliableSites/D/";

        String resultDirS = "/Volumes/VMap2_Fei/reliableSites/ABD_intersect";

        this.mkReliable(abDirS, abdDirS, resultDirS);
        this.mkReliable(dDirS, abdDirS, resultDirS);

    }

    private void mkReliable (String source1DirS, String source2DirS, String resultDirS) {
        List<File> fList = IOUtils.getFileListInDirEndsWith(source1DirS, ".gz");
        Collections.sort(fList);
        fList.parallelStream().forEach(f -> {
            String source2FileS = new File (source2DirS, f.getName().split("_")[0]+"_ABD_reliable.txt.gz").getAbsolutePath();
            String outputFileS = new File (resultDirS, f.getName().split("_")[0]+"_intersect_reliable.txt.gz").getAbsolutePath();
            try {
                BufferedReader br1 = IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedReader br2 = IOUtils.getTextGzipReader(source2FileS);
                BufferedWriter bw = IOUtils.getTextGzipWriter(outputFileS);
                String temp1 = br1.readLine();
                String temp2 = br2.readLine();
                bw.write(temp1);
                bw.newLine();
                int cnt = 0;
                while ((temp1 = br1.readLine()) != null) {
                    temp2 = br2.readLine();
                    bw.write(String.valueOf(Integer.parseInt(temp1)*Integer.parseInt(temp2)));
                    bw.newLine();
                    cnt++;
                    if (cnt%10000000 == 0) System.out.println(cnt);
                }
                bw.flush();
                bw.close();
                br1.close();
                br2.close();
                System.out.println(f.getName());
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void mkReliableGenotypeSite () {
        String abInDirS = "/Volumes/VMap2_Fei/popdep_vmap2/AB/";
        String abdInDirS = "/Volumes/VMap2_Fei/popdep_vmap2/ABD/";
        String dInDirS = "/Volumes/VMap2_Fei/popdep_vmap2/D/";

        String abOutDirS = "/Volumes/VMap2_Fei/reliableSites/AB";
        String abdOutDirS = "/Volumes/VMap2_Fei/reliableSites/ABD";
        String dOutDirS = "/Volumes/VMap2_Fei/reliableSites/D";

        String abSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abPopDep_sample.txt";
        String abdSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abdPopDep_sample.txt";
        String dSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/dPopDep_sample.txt";

        double depthStart = 2;
        double depthEnd = 8;
        double SDStart = 2;
        double SDEnd = 8;
        int nBin = 100;
        double proportionOfSite = 0.7;

        this.mkReliable(abInDirS, abOutDirS, abSampleFileS, depthStart, depthEnd, SDStart, SDEnd, nBin, proportionOfSite);

        depthStart = 3;
        depthEnd = 15;
        SDStart = 3;
        SDEnd = 8;
        this.mkReliable(abdInDirS, abdOutDirS, abdSampleFileS, depthStart, depthEnd, SDStart, SDEnd, nBin, proportionOfSite);

        depthStart = 3;
        depthEnd = 17;
        SDStart = 3;
        SDEnd = 10;
        this.mkReliable(dInDirS, dOutDirS, dSampleFileS, depthStart, depthEnd, SDStart, SDEnd, nBin, proportionOfSite);

    }

    private void mkReliable (String inDirS, String outDirS, String sampleFileS, double depthStart, double depthEnd, double SDStart, double SDEnd, int nBin, double proportionOfSite) {
        RowTable<String> t = new RowTable<>(sampleFileS);
        Grid gr = new Grid(depthStart, depthEnd, SDStart, SDEnd, nBin);
        for (int i = 0; i < t.getRowNumber(); i++) {
            gr.addXY(t.getCellAsDouble(i,1), t.getCellAsDouble(i,2));
        }
        gr.buildHashMap();
        int indexThresh = gr.getOrderIndexOfProportionOfSite(proportionOfSite);
        List<File> fList = IOUtils.getFileListInDirEndsWith(inDirS,".gz");
        fList.stream().forEach(f -> {
            String outfileS = f.getName().replaceFirst("popdep_vmap2.txt.gz", "reliable.txt.gz");
            outfileS = new File(outDirS, outfileS).getAbsolutePath();
            try {
                BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                bw.write("IfReliable(0/1)");
                bw.newLine();
                List<String> l = new ArrayList<>();
                String temp = br.readLine();
                double x;
                double y;
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    x = Double.parseDouble(l.get(1));
                    y = Double.parseDouble(l.get(2));
                    if (gr.isHighDensity(x, y, indexThresh)) {
                        bw.write("1");
                    }
                    else {
                        bw.write("0");
                    }
                    bw.newLine();
                    cnt++;
                    if (cnt%10000000 == 0) System.out.println(cnt);
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName());
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }



    public void densityFilter () {
        String abSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abPopDep_sample.txt";
        String abdSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/abdPopDep_sample.txt";
        String dSampleFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/dPopDep_sample.txt";
        double proportionOfSite = 0.70;
        this.densityFilterPlotAB(proportionOfSite, abSampleFileS);
        this.densityFilterPlotABD(proportionOfSite, abdSampleFileS);
        this.densityFilterPlotD(proportionOfSite, dSampleFileS);
    }



    public void densityFilterPlotD(double proportionOfSite, String inputFileS) {
        String beforeFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/beforeAfter/d_before_scatter.pdf";
        String afterFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/beforeAfter/d_after_scatter.pdf";
        RowTable<String> t = new RowTable<>(inputFileS);
        double depthStart = 3;
        double depthEnd = 17;
        double SDStart = 3;
        double SDEnd = 10;
        int nBin = 100;
        Grid gr = new Grid(depthStart, depthEnd, SDStart, SDEnd, nBin);
        for (int i = 0; i < t.getRowNumber(); i++) {
            gr.addXY(t.getCellAsDouble(i,1), t.getCellAsDouble(i,2));
        }
        gr.buildHashMap();
        TDoubleArrayList xList = new TDoubleArrayList();
        TDoubleArrayList yList = new TDoubleArrayList();
        int indexThresh = gr.getOrderIndexOfProportionOfSite(proportionOfSite);
        for (int i = 0; i < t.getRowNumber(); i++) {
            double x = t.getCellAsDouble(i, 1);
            double y = t.getCellAsDouble(i, 2);
            if (!gr.isHighDensity(x, y , indexThresh)) continue;
            xList.add(x);
            yList.add(y);
        }

        ScatterPlot s = new ScatterPlot(t.getColumnAsDoubleArray(1), t.getColumnAsDoubleArray(2));
        s.setTitle("D_before");
        s.setColor(255,0, 0, 5);
        s.setXLim(0, 20);
        s.setYLim(0, 12);
        s.setXLab("Mean of depth");
        s.setYLab("SD of depth");
        s.saveGraph(beforeFileS);
        s = new ScatterPlot(xList.toArray(), yList.toArray());
        s.setTitle("D_after");
        s.setColor(255,0, 0, 5);
        s.setXLim(0, 20);
        s.setYLim(0, 12);
        s.setXLab("Mean of depth");
        s.setYLab("SD of depth");
        s.saveGraph(afterFileS);
    }
    public void densityFilterPlotABD(double proportionOfSite, String inputFileS) {
        String beforeFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/beforeAfter/abd_before_scatter.pdf";
        String afterFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/beforeAfter/abd_after_scatter.pdf";
        RowTable<String> t = new RowTable<>(inputFileS);
        double depthStart = 3;
        double depthEnd = 15;
        double SDStart = 3;
        double SDEnd = 8;
        int nBin = 100;
        Grid gr = new Grid(depthStart, depthEnd, SDStart, SDEnd, nBin);
        for (int i = 0; i < t.getRowNumber(); i++) {
            gr.addXY(t.getCellAsDouble(i,1), t.getCellAsDouble(i,2));
        }
        gr.buildHashMap();
        TDoubleArrayList xList = new TDoubleArrayList();
        TDoubleArrayList yList = new TDoubleArrayList();
        int indexThresh = gr.getOrderIndexOfProportionOfSite(proportionOfSite);
        for (int i = 0; i < t.getRowNumber(); i++) {
            double x = t.getCellAsDouble(i, 1);
            double y = t.getCellAsDouble(i, 2);
            if (!gr.isHighDensity(x, y , indexThresh)) continue;
            xList.add(x);
            yList.add(y);
        }

        ScatterPlot s = new ScatterPlot(t.getColumnAsDoubleArray(1), t.getColumnAsDoubleArray(2));
        s.setTitle("ABD_before");
        s.setColor(255,0, 0, 5);
        s.setXLim(0, 20);
        s.setYLim(0, 12);
        s.setXLab("Mean of depth");
        s.setYLab("SD of depth");
        s.saveGraph(beforeFileS);
        s = new ScatterPlot(xList.toArray(), yList.toArray());
        s.setTitle("ABD_after");
        s.setColor(255,0, 0, 5);
        s.setXLim(0, 20);
        s.setYLim(0, 12);
        s.setXLab("Mean of depth");
        s.saveGraph(afterFileS);
    }

    public void densityFilterPlotAB(double proportionOfSite, String inputFileS) {
        String beforeFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/beforeAfter/ab_before_scatter.pdf";
        String afterFileS = "/Users/feilu/Documents/analysisL/production/vmap2/depth/plot/beforeAfter/ab_after_scatter.pdf";
        RowTable<String> t = new RowTable<>(inputFileS);
        double depthStart = 2;
        double depthEnd = 8;
        double SDStart = 2;
        double SDEnd = 8;
        int nBin = 100;
        Grid gr = new Grid(depthStart, depthEnd, SDStart, SDEnd, nBin);
        for (int i = 0; i < t.getRowNumber(); i++) {
            gr.addXY(t.getCellAsDouble(i,1), t.getCellAsDouble(i,2));
        }
        gr.buildHashMap();
        TDoubleArrayList xList = new TDoubleArrayList();
        TDoubleArrayList yList = new TDoubleArrayList();
        int indexThresh = gr.getOrderIndexOfProportionOfSite(proportionOfSite);
        for (int i = 0; i < t.getRowNumber(); i++) {
            double x = t.getCellAsDouble(i, 1);
            double y = t.getCellAsDouble(i, 2);
            if (!gr.isHighDensity(x, y , indexThresh)) continue;
            xList.add(x);
            yList.add(y);
        }
        ScatterPlot s = new ScatterPlot(t.getColumnAsDoubleArray(1), t.getColumnAsDoubleArray(2));
        s.setTitle("AB_before");
        s.setColor(255,0, 0, 5);
        s.setXLim(0, 8);
        s.setYLim(0, 8);
        s.setXLab("Mean of depth");
        s.setYLab("SD of depth");
        s.saveGraph(beforeFileS);
        s = new ScatterPlot(xList.toArray(), yList.toArray());
        s.setTitle("AB_after");
        s.setColor(255,0, 0, 5);
        s.setXLim(0, 8);
        s.setYLim(0, 8);
        s.setXLab("Mean of depth");
        s.setYLab("SD of depth");
        s.saveGraph(afterFileS);
    }

    /**
     * By aoyue
     */
    public void findPopDepMode () {
        double abDepthMode = 5.144204;
        double abdDepthMode = 9.607552;
        double dDepthMode = 11.91677;
        double abSDMode = 4.226852;
        double abdSDMode = 5.162204;
        double dSDMode = 7.056306;
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
