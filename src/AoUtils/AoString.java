package AoUtils;

import java.util.List;

public class AoString {
    public AoString(){

    }
//String name = new File(outPath).getName();
    //String chr = name.substring(name.indexOf("chr")+3,name.indexOf("chr")+5);
    //String pop1pop2 = name.substring(0,name.indexOf("_chr"));

    public static String listToString(List list, String separator) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < list.size(); i++) {
            if (i == list.size() - 1) {
                sb.append(list.get(i));
            } else {
                sb.append(list.get(i));
                sb.append(separator);
            }
        }
        return sb.toString();
    }
}
