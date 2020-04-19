package daxing.filterSNP;

import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class DepthInfocp {

    private List<Dot> dotList;

    public DepthInfocp(String inputFile){
        this.initialize(inputFile);
    }

    private void initialize(String inputFile){
        List<Dot> dotList=new ArrayList<>();
        try(BufferedReader br= IOUtils.getTextGzipReader(inputFile)){
            String line;
            List<String> lineList;
            short chr;
            int pos;
            double depth;
            double sd;
            br.readLine();
            while ((line=br.readLine())!=null){
                lineList= PStringUtils.fastSplit(line);
                chr=Short.parseShort(lineList.get(0));
                pos=Integer.parseInt(lineList.get(1));
                depth=Double.parseDouble(lineList.get(2));
                sd=Double.parseDouble(lineList.get(3));
                dotList.add(new Dot(chr, pos, depth, sd));
            }
            Collections.sort(dotList);
            this.dotList=dotList;
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public List<Dot> getDotList() {
        return dotList;
    }
}
