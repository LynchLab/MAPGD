The row class servers as a cheep and easy way to pass information between different routines within mapgd. It is intended to be easy for new developers to implement, to help ensure that new mapgd commands can benefit from multi-threaded and multi-node execution, and to serve as a convenient way for users and developers alike to interact with the output of mapgd commands.

empty	| KEY NAME      | KEY NAME      | KEY NAME      
	| key info.	| key info.	| key info.	
------- | ------------- | ------------- | ---------------
row 	|   datum	|    datum	|   datum	
row 	|   datum	|    datum	|   datum	
row 	|   datum	|    datum	|   datum	

First Header  | Second Header
------------- | -------------
Content Cell  | Content Cell
Content Cell  | Content Cell

For instance, the output of a text command can look like this:

```			
CHROM  |POS    |LOCUS                                         
(id#)  |(bp)   |(An array of quartets!)
       |       |IND\_1         IND\_2         IND\_3         IND\_4
--------------------------------------------------------------------
0       1       0/0/0/0        0/0/0/1        199/0/0/0      2/0/0/1
0       2       0/1/0/0        0/1/0/0        0/0/0/0        2/0/0/1
```
				
 new row can b

```C++
void row.get(key, void* dst);
void row.put(key, void* src);

void * row.get(key get_key);
void * row.put(key put_key);

void row.read(std::istream in);
void row.write(std::ostream out);

key row.get_keys();
void row.set_keys(container <key> keys);
void row.set_keys(std::istream in);
```

```C++
keyid_t key (std::istream in);				
keyid_t key (char* name, std::string desc);		
keyid_t key (char* name, std::string desc, format); 	
```

```C++
char* key.get_name();
keyid_t key.get_id();
std::string key.get_desc();
```

```C++
format(data);

void format.get_key(void* src, char* dst);
void format.put_key(void* src, char* dst);
void format.get_value(void* src, char* dst);
void format.put_value(char* src, void* dst);

void format.get_value_bin(void* src, char* dst);
void format.put_value_bin(char* src, void* dst);
void format.get_value_txt(void* src, char* dst);
void format.put_value_txt(char* src, void* dst);
```
