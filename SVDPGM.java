package svdpgm;
import java.io. *;

import svdpgm.PGMIO;
import Jama.Matrix;
import Jama.SingularValueDecomposition;

public class SVDPGM {

	/**
	 * Writes a matrix to file delimiting by whitespace
	 * 
	 * @param matrix		Matrix to be written to file
	 * @param filename		Filename of file
	 */
	public static void writeMatrixToFile(Matrix matrix, String filename){
		try (Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename), "utf-8"))) {
			for (int i = 0; i < matrix.getRowDimension(); i++){
				for (int j = 0; j < matrix.getColumnDimension(); j++){
					writer.write(String.format("%f", matrix.get(i, j)));
					if( j != matrix.getColumnDimension()-1){
						writer.write(" ");
					}
				}
				writer.write("\n");
			}
			writer.close();
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Generates a matrix corresponding to some given image
	 * 
	 * @param filename		filename of image to be converted to matrix
	 * @return
	 * @throws FileNotFoundException 
	 * @throws IOException
	 */
	public static Matrix generateMatrix(String filename) throws IOException{
		int[][] imageIntArray;
		double [][] imageDoubleArray;
		Matrix imageMatrix;
		File image = new File(filename);

		imageIntArray = PGMIO.read(image);
		imageDoubleArray = convert2DIntTo2DDouble(imageIntArray);
		imageMatrix = new Matrix(imageDoubleArray);
		
		return imageMatrix;
	}
	
	/**
	 * Generates images Pseudoinverse.pgm and IHat.pgm, corresponding to the pseudoinverse of a given image and the result of multiplying the original with the psuedoinverse
	 * 
	 * @param filename		image to generate pseudoinverse for
	 * @throws IOException
	 * @return 				matrix corresponding to pseudoinverse of image
	 */
	public static Matrix generatePseudoinverse(String filename) throws IOException{
		Matrix imageMatrix = generateMatrix(filename);
		
		return generatePseudoinverse(imageMatrix);

	}
	
	/**
	 * Generates images Pseudoinverse.pgm and IHat.pgm, corresponding to the pseudoinverse of the image represented by a matrix 
	 * and the result of multiplying the original with the psuedoinverse
	 * 
	 * @param matrix		matrix to generate pseudoinverse for
	 * @throws IOException
	 * @return 				matrix corresponding to pseudoinverse of matrix
	 */
	public static Matrix generatePseudoinverse(Matrix matrix) throws IOException{

		SingularValueDecomposition SVD = new SingularValueDecomposition(matrix);
		Matrix U = SVD.getU();
		Matrix UTranspose = U.transpose();
		
		Matrix D = SVD.getS();
		double [][] DArray = D.getArray();
		double [][] DPlusArray = new double [DArray.length][];
		
		for(int i = 0; i < DArray.length; i++){
			DPlusArray[i] = new double [DArray[i].length];
			for (int j = 0; j < DArray.length; j++){
				if(DArray[i][j] > 0){
					DPlusArray[i][j] = 1/DArray[i][j];
				}
				else{
					DPlusArray[i][j] = 0;
				}
			}
		}
		Matrix DPlus = new Matrix(DPlusArray);
		
		Matrix V = SVD.getV();
		
		Matrix pseudoinverse = V.times(DPlus).times(UTranspose);
		generateImage(pseudoinverse, "Pseudoinverse.pgm");;
		
		return pseudoinverse;

	}
	
	/**
	 * Generates image A"p".pgm, where "p" is replaced by the value of p, as an approximation of a given image via SVD
	 * and returns the matrix corresponding to the approximation.
	 * 
	 * @param filename		filename of image to be approximated
	 * @param p				number of singular values to use in generating approximation, at least 1
	 * @throws IOException
	 * @return				matrix corresponding to image approximation
	 */
	public static Matrix generateApproximation(String filename, int p) throws IOException{
		Matrix imageMatrix = generateMatrix(filename);
		
		return generateApproximation(imageMatrix, p);
	}
	
	/**
	 * Generates image A"p".pgm, where "p" is replaced by the value of p, as an approximation of a given image via SVD
	 * and returns the matrix corresponding to the approximation.
	 * 
	 * @param imageMatrix	matrix generated from image
	 * @param p				number of singular values to use in generating approximation, at least 1
	 * @throws IOException
	 * @return				matrix corresponding to image approximation
	 */
	public static Matrix generateApproximation(Matrix imageMatrix, int p) throws IOException{

		SingularValueDecomposition SVD = new SingularValueDecomposition(imageMatrix);
		Matrix U = SVD.getU();
		Matrix UTranspose = U.transpose();
		double[][] UTransposeArray = UTranspose.getArray();
		double[] S = SVD.getSingularValues();
		Matrix V = SVD.getV();
		Matrix VTranspose = V.transpose();
		double[][] VArray = VTranspose.getArray();
		
		double [][] u1Array = new double [1][1];
		u1Array[0] = UTransposeArray[0];
		Matrix u1 = new Matrix(u1Array);
		u1 = u1.transpose();
				
		double [][] v1Array = new double [1][1];
		v1Array[0] = VArray[0];
		Matrix v1 = new Matrix(v1Array);
		
		Matrix A = u1.times(v1).times(S[0]);
		
		for(int k = 1; k < p; k++){
			double [][] uArray = new double [1][1];
			uArray[0] = UTransposeArray[k];
			Matrix u = new Matrix(uArray);
			u = u.transpose();
					
			double [][] vArray = new double [1][1];
			vArray[0] = VArray[k];
			Matrix v = new Matrix(vArray);
			
			A = A.plus(u.times(v).times(S[k]));
		}
		
		generateImage(A, "A" + Integer.toString(p) +".pgm");
		return A;
	}
	
	/**
	 * Generates images U.pgm, D.pgm, VTranspose.pgm from SVD of given file
	 * 
	 * @param filename	name of file to generate SVD for.
	 * @throws IOException
	 */
	public static void generateSVD(String filename) throws IOException{
		Matrix imageMatrix = generateMatrix(filename);
		
		SingularValueDecomposition SVD = new SingularValueDecomposition(imageMatrix);
		Matrix U = SVD.getU();
		Matrix D = SVD.getS();
		Matrix V = SVD.getV();
		V = V.transpose();
		
		generateImage(U, "U.pgm");
		generateImage(D, "D.pgm");
		generateImage(V, "VTranspose.pgm");
	}
	
	/**
	 * Generates an image from a matrix
	 * 
	 * @param imageMatrix	Matrix of image info
	 * @param filename		Name for image
	 * @throws IOException
	 */
	public static void generateImage(Matrix imageMatrix, String filename) throws IOException{
		double[][] imageDoubleArray = imageMatrix.getArray();
		int[][] imageIntArray = mapArray(imageDoubleArray);
		File file = new File(filename);
		PGMIO.write(imageIntArray, file);
	}
	
	 /**
	  * Maps linearly the elements from a 2D double array to range [0, 255]
	  * 
	  * @param array	the array to be mapped
	  * @return			the array that has had elements mapped
	  */
	public static int[][] mapArray(double[][] array){
		double max = array[0][0];
		double min = array[0][0];
		for (int i = 0; i < array.length; i++){
			for (int j = 0; j < array[i].length; j++){
				if(array[i][j] > max){
					max = array[i][j];
				}
				if(array[i][j] < min){
					min = array[i][j];
				}
			}
		}
		
		int[][] intArray = new int[array.length][];
		for (int i = 0; i < array.length; i++){
			intArray[i] = new int[array[i].length];
			for (int j = 0; j < array[i].length; j++){
				intArray[i][j] = (int)Math.round(255*((array[i][j] - min) / (max - min)));
			}
		}
		return intArray;
	}
	
	/**
	 * Converts a 2D int array into a 2D float array
	 * 
	 * @param int2DArray
	 * @return
	 */
	public static double[][] convert2DIntTo2DDouble(int[][] int2DArray){
		double[][] double2DArray = new double[int2DArray.length][];
		for (int i = 0; i < int2DArray.length; i++){
			double2DArray[i] = new double[int2DArray[i].length];
			for (int j = 0; j < int2DArray[i].length; j++){
				double2DArray[i][j] = int2DArray[i][j];
			}
		}
		return double2DArray;
	}

}
