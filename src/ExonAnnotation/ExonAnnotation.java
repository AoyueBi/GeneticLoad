package ExonAnnotation;

import AoUtils.AoFile;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class ExonAnnotation {

    String exonAnnotationFileS = "";

    //************** 建立对应的 HashMap *********************//
    private static HashMap<String, String> hmIDIfRefisAnc = null;
    private static HashMap<String, String> hmIDLoadGroup = null;
    private static HashMap<String, String> hmTaxaGenomeType = null;


    public ExonAnnotation(){
        this.readExonAnnotation();
    }

    public void readExonAnnotation(){



        try {
            BufferedReader br = AoFile.readFile(exonAnnotationFileS);
            String header = br.readLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
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
