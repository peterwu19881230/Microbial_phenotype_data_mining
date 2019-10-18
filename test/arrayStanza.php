<?php
//Goal: Write a function that is able to parse all the element in an array or associative array into one string delimited by space

//parse the .CSV file into an array
$ECKarray = array_map('str_getcsv', file('ECK.csv'));

//These 2 lines call the function array_stanza($array) and ouput the result
$result=array_stanza($ECKarray);
var_dump($result);


//Recursive function. Will re-enter the function if the next level is still an array. If it is string, they will be concatnated into $stanza and output.
function array_stanza($array){

	
	
	if(!is_string($array[0])){
		
		//Must declare $nextArray before putting it in array_merge  			
		$nextArray=array();
		
		foreach ($array as $new_array){
		/*
		Warning: Be aware that array_merge() overwrites the non-numeric keys. 
		*/
		$nextArray=array_merge($nextArray,$new_array);
		
		}
		array_stanza($nextArray);		
		
	}
	else{
		foreach ($array as $string){
		$stanza=$stanza." ".$string; 		
		}
		
	}
	/*
	echo($stanza); works and I tried to only put several characters for the return value, but the return value always gets me a NULL.	
	*/
	return $stanza;
	
	

	
	
}






