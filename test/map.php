<?php
/*
Goal: Map ECK strains with JW number and find those that don't have JW number. 

About the code: ECK.csv and JW.csv containing strain names were used. 2 Foreach loops were used: 1st is to to map the JW number for ECK strains in Nichols' paper, 
2nd is to find ECK strains in Nichols' paper that don't have JW number. The result is in JWInNicholsPaper.txt for further usage.
Statistically, 88 strains were found not having JW number. But it's wierd because the paper claims that more than 200 don't have ECK number. 
*/

/* This function is an old way I found in php website. I now use just array_map()
function capture_matrix($filename){

	$row = 1;
	if (($handle = fopen("$filename", "r")) !== FALSE) {
   	 while (($data = fgetcsv($handle, 1000, ",")) !== FALSE) {
    	    $num = count($data);
        	$row++;
        	for ($c=0; $c < $num; $c++) {                       
         	$array[$row][$c]=$data[$c];   
        	}
   
    	}
    	fclose($handle);
    
    	return $array;
}

}
*/

//parse the .CSV file into an array
$ECKarray = array_map('str_getcsv', file('ECK.csv'));
$JWarray = array_map('str_getcsv', file('JW.csv'));




//$ECKarray=capture_matrix("ECK.csv");
//$JWarray=capture_matrix("JW.csv");

//I reluctantly have to create a new array for ECK.csv because the original captured array was an associative array (each strain is actually a subarray with the key being "0")
foreach ($ECKarray as $ECK)
{
$NewECKarray[]=$ECK[0];
}





//Compare the 2 arrays:
$n=0; $m=0;

//This line is for creating an empty file (will overwrite the old one)
file_put_contents("JWInNicholsPaper.txt","");

foreach ($JWarray as $JW)
{	
		
		
		
		
		
		
		$JW[0]=preg_quote($JW[0]);
		if(count(preg_grep("($JW[0])", $NewECKarray))!=0)
		{
		file_put_contents("JWInNicholsPaper.txt","$JW[0] $JW[2]\n",FILE_APPEND);
		$n++;		
		}				
}
file_put_contents("JWInNicholsPaper.txt","$n strains were found\n",FILE_APPEND);


//Compare 2 arrays. Reset $n to 0.
$n=0;

//I reluctantly have to create a new array for ECK.csv because the original captured array was an associative array (each strain is actually a subarray with the key being "0")
foreach ($JWarray as $JW)
{
$JWarrayFirstColumn[]=$JW[0];
}

foreach ($NewECKarray as $ECK)
{	
	//Delete the expression of the deleted gene after ECK number in $ECKarray
	$modifiedECK=preg_replace("(-.*)",'',$ECK);
	
	
	
		if(count(preg_grep("($modifiedECK)", $JWarrayFirstColumn))==0)
		{
		file_put_contents("JWInNicholsPaper.txt","$modifiedECK no JW!\n",FILE_APPEND);
		echo $modifiedECK." was not found to have JW number\n";
		$n++;
		}
}
echo "->$n strains were not found to have JW number\n";


