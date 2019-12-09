/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package AoUtils;

import utils.IOUtils;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author Aoyue
 */
public class AoMath {

    public AoMath() {

    }

    /**
     * 过滤DAF_ABD等于0或者1的位点
     *
     */
    public void filterValue(String infileDirS, String outfileDirS) {

//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/004_merge";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/104_feiResult/005_filterDAF_ABD";
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach((File f) -> {

            try {
                String infileS = f.getAbsolutePath();
                String outfileS = new File(outfileDirS, f.getName().replaceFirst(".txt", "_filterGERPandPhyloP.txt")).getAbsolutePath();

                BufferedReader br = null;
                if (infileS.endsWith(".txt")) {
                    br = IOUtils.getTextReader(infileS);
                } else if (infileS.endsWith(".txt.gz")) {
                    br = IOUtils.getTextGzipReader(infileS);
                }
                BufferedWriter bw = null; // IOUtils.getTextGzipWriter(outfileS);
                if (outfileS.endsWith(".txt")) {
                    bw = IOUtils.getTextWriter(outfileS);
                } else if (outfileS.endsWith(".txt.gz")) {
                    bw = IOUtils.getTextGzipWriter(outfileS);
                }
                String temp = null;
                String header = br.readLine(); //读表头
                bw.write(header);
                bw.newLine();
                int cnt = 0;
                List<String> l = new ArrayList();
                String goalValue1 = null;
                String goalValue2 = null;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    StringBuilder sb = new StringBuilder();
                    cnt++;
//ID	Chr	Pos	Ref	Alt	Major	Minor	Maf	AAF_ABD	AAF_AB	Transcript	Region	Variant_type	SIFT_score	Ancestral	DAF	DAF_ABD	DAF_AB	Gerp	PhyloP
                    goalValue1 = l.get(18); //此处需要修改，目标值
                    goalValue2 = l.get(19);
                    if (goalValue1.startsWith("N")) {
                        continue;
                    }
                    if (goalValue2.startsWith("N")) {
                        continue;
                    }

                    double value1 = Double.parseDouble(goalValue1);
                    double value2 = Double.parseDouble(goalValue2);
//                    if (value1 == 0 || (value1 == 1)) {
//                        continue;
//                    }
                    if (value1 < 1 || (value2 < 0.5)) {
                        continue;
                    }
                    bw.write(temp);
                    bw.newLine();

                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getAbsolutePath() + " is completed");
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

        });
        /*==================================== 测试用 =============================================*/

    }

}
