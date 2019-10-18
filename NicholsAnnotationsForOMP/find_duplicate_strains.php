<?php
//Goal: Find all the strains in Nichols' paper that appeared more than twice (Which shouldn't )



#this opens the strain mmc file
		$filename='ECK.csv';
		echo"opening $filename\n";
		$fh=fopen($filename, "r");

//parse the .CSV file into an array
$ExcelArray = array_map('str_getcsv', file('ECK.csv'));

//This foreach is to make $ExcelArray into another simplied array, instead of having subarray in it for each string.
foreach($ExcelArray as $SubArray){

$ExcelSimplifiedArray[]=$SubArray[0];

}


$UniqueExcelSimplifiedArray=array_unique($ExcelSimplifiedArray);


/*
array_count_values analyzes the times each element appears in an array. The result is an array.
http://stackoverflow.com/questions/13633954/how-do-i-count-occurrence-of-duplicate-items-in-array
*/
$result=array_count_values($ExcelSimplifiedArray);

$n=0;
foreach ($result as $strain=>$TimesUsed){

	If($TimesUsed>=2){
	echo $strain." was used ".$TimesUsed." times\n";
	$n++;
	}


}
echo $n." strains were used at least twice.\n";