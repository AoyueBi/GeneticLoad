/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import pgl.infra.position.ChrPos;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 *
 * @author Aoyue
 */
public class Mode {

    /*==================================== 测试用 =============================================*/
    //*********************************** START1 ***********************************//
    //***************************** step one : 确定其倍性，根据倍性计算 ****************************//

    public Mode() {

    }

    public void mergeTxt(){
        String infileDirS = "";
        String outfileS = "";
        AoFile.mergeTxt(infileDirS,outfileS);
    }

    public void runParallele_listFile(){ //本地运行常用
        String infileDirS = "/Users/Aoyue/Documents/out";
        String outfileDirS = "/Users/Aoyue/Documents/out1";
        List<File> fsList = IOUtils.getVisibleFileListInDir(infileDirS);
        fsList.stream().forEach(f -> {
            String infileS = f.getAbsolutePath();
            String outfileS = null;
            if (infileS.endsWith(".txt")) {
                outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "_RH_2Mwindow_1Mstep.txt.gz").getAbsolutePath();
            } else if (infileS.endsWith(".txt.gz")) {
                outfileS = new File(outfileDirS, f.getName().split(".txt.gz")[0] + "_RH_2Mwindow_1Mstep.txt.gz").getAbsolutePath();
            }
            new Bin().calwindowstep_ResidualHeterozygosity(infileS,2000000,1000000,outfileS);
            System.out.println(f.getName() + "\tis completed at " + outfileS);
        });
    }



    public void runJarParallele(){ //一次性将所有的jar都运行上
        String infileDirS = "";
        String outfileDirS ="";
        String logDirS = "";
//        String[] chrArr = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
//        String[] chrArr = {"1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"};
//        String[] chrArr = {"1D","2D", "3D", "4D", "5D", "6D","7D"};
//        String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B", "1D", "2D", "3D",  "4D", "5D", "6D", "7D"};
        String[] chrArr = {"1A", "2A", "3A","4A", "5A", "6A", "7A", "1B", "2B", "3B", "4B", "5B", "6B", "7B"};

//        String[] chrArr = {"001","002","003","004","005","006","007","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022","023","024","025","026","027","028","029","030","031","032","033","034","035","036","037","038","039","040","041","042"};
//        String[] chrArr ={"001","002","003","004","007","008","009","010","013","014","015","016","019","020","021","022","025","026","027","028","031","032","033","034","037","038","039","040"};
//        String[] chrArr ={"005","006","011","012","017","018","023","024","029","030","035","036","041","042"};

        for (int i = 0; i < chrArr.length; i++) {
            String infileS = new File(infileDirS,"chr" + chrArr[i] + "_vmap2.1_heter_SNPbased_Cultivar.vcf.gz").getAbsolutePath();
            String outfileS = new File(outfileDirS,"chr" + chrArr[i] + "_vmap2.1_heter_SNPbased_Cultivar_chrposGenotype.txt.gz").getAbsolutePath();
            String logfileS = new File(logDirS,"log_" + new File(outfileS).getName().split(".gz")[0]).getAbsolutePath(); //不管是不是gz结尾，我们只取gz前的部分，妙！
            System.out.println("nohup java -jar 034_mkindividualVCFtoChrPosGenotype.jar " + infileS + " " + outfileS + " > " + logfileS  + " 2>&1 &" );
        }

    }


    public void vcfParallel() {
        double extractRatio = 0;
        String infileDirS = "";
        String outfileDirS = "";

        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        Collections.sort(fsList);
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = null;
                BufferedReader br = null;
                if (infileS.endsWith(".vcf")) {
                    br = IOUtils.getTextReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().split(".vcf")[0] + "_subset.vcf.gz").getAbsolutePath();
                } else if (infileS.endsWith(".vcf.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                    outfileS = new File(outfileDirS, f.getName().split(".vcf.gz")[0] + "_subset.vcf.gz").getAbsolutePath();
                }
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                String temp = null;
                int cnttotal = 0;
                int cntsubset = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp);
                        bw.newLine();
                    } else {
                        cnttotal++;
                        double r = Math.random();
                        double ratio = extractRatio;
                        if (r > ratio) {
                            continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                        }
                        List<String> l = PStringUtils.fastSplit(temp);
                        if (l.get(3).contains(",")) {
                            continue; // 第3列是alt的信息，若有2个等位基因，则去除这一行
                        }
                        bw.write(temp);
                        bw.newLine();
                        cntsubset++;
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName() + "\twith " + cnttotal + " bp has a subset of\t" + cntsubset + "\tis completed at " + outfileS);
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    //本模板使用多线程流
    public void txtParallel() {
        String infileDirS = "";
        String outfileDirS = "";
        List<File> fsList = AoFile.getFileListInDir(infileDirS);
        fsList.parallelStream().forEach(f -> {
            try {
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS, f.getName().split(".txt")[0] + "_subset.txt.gz").getAbsolutePath();
                BufferedReader br = AoFile.readFile(infileS);
                BufferedWriter bw = AoFile.writeFile(outfileS);
                String header = br.readLine();
                String temp = null;
                List<String> l = new ArrayList<>();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);

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

    public void vcfSinglethread() {


        try {
            String infileS = "";
            String outfileS = "";
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            int cnt = 0;
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                StringBuilder sb = new StringBuilder();
                if (temp.startsWith("#")) {
                    continue;
                }
                cnt++;
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }


        List<ChrPos> l = new ArrayList();
        String chr = null;
        String pos = null;
        l.add(new ChrPos(Short.valueOf(chr), Integer.valueOf(pos)));

    }


    /**
     *  ################################### 初
     */

    public void txtSinglethread(){
        String infileS = "";
        String outfileS = "";
        try {
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String header = br.readLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;

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

    public void TryCatch(){
        try {
            String infileS = "";
            String outfileS = "";
            BufferedReader br = AoFile.readFile(infileS);
            BufferedWriter bw = AoFile.writeFile(outfileS);
            String temp = null;
            String header = br.readLine();
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                cnt++;

            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        try{

        }catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    public void testifD() {
        boolean ifd = false;
        String chr = "chr003.vmap2...".substring(3, 6);
        //根据染色体号进行AB还是D的判断
        String[] db = {"5", "6", "11", "12", "17", "18", "23", "24", "29", "30", "35", "36", "41", "42"};
        Arrays.sort(db);
        if (Arrays.binarySearch(db, chr) > -1) { //说明是属于D的
            ifd = true;
        }
        try {
            String infileS = "";
            String outfileS = "";
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;

            if (ifd == false) {
                bw.write(temp + "\tAncestral\tDaf\tDaf_ABD\tDaf_AB");
                bw.newLine();
            } else if (ifd == true) {
                bw.write(temp + "\tAncestral\tDaf\tDaf_ABD\tDaf_D");
                bw.newLine();
            }

            while ((temp = br.readLine()) != null) {
                cnt++;

            }
            br.close();
            System.out.println();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

}
