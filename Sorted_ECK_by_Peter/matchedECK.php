<?php
//Goal: parse Keio_strains.csv into an associative array 


$JWarray = array_map('str_getcsv', file('Keio_strains_with_verified_JW.csv'));

var_dump($JWarray);
