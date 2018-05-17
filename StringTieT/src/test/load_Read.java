package test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.Vector;

public class load_Read {
	public static Vector<String> read_vector=new Vector<String>();
	/**
	 * 从输入文件中解析出reads信息，并将reads信息保存在read_vector中
	 * @throws IOException
	 */
	   public static  void load_reads() throws IOException{
		   File file=new File("sim130bp.fa");
		   BufferedReader reader=null;
		  // read_vector;
		   try {
			reader=new BufferedReader(new FileReader(file));
			String tempString=null;
			while((tempString=reader.readLine())!=null){
				tempString=reader.readLine();
				read_vector.addElement(tempString);
				//System.out.println(tempString);
			}
//			Iterator iterator=read_vector.iterator();
//			while(iterator.hasNext()){
//				System.out.println(iterator.next());
//			}
//			System.out.println(length);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	   }


}