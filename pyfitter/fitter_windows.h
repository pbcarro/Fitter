//Program setup functions
__declspec(dllexport) int Initialize_Stuff (double **/*ETArray*/, int */*CatTransitions*/, int */*DictTransitions*/, double */*FileDelta*/, int */*StatePoints*/, struct Level **/*DictionaryIn*/, struct Transition **/*CatalogIn*/);
__declspec(dllexport) int Load_ETau_File (char */*FileName*/, double **/*X*/, double */*FileDelta*/, int */*StatePoints*/, int */*StateCount*/);
__declspec(dllexport) int Load_Base_Catalog (char */*FileName*/, struct Transition **/*BaseCatalog*/,  int /*Verbose*/);
__declspec(dllexport) int Load_Base_Catalog_Dictionary (char */*FileName*/, struct Level **/*DictIn*/,  int /*Verbose*/);
__declspec(dllexport) int Load_Exp_File  (char */*FileName*/, double **/*X*/, double **/*Y*/, int /*Verbose*/);

//Frequency predicting functions
__declspec(dllexport) double Get_Kappa (double /*A*/, double /*B*/, double /*C*/);  
__declspec(dllexport) int Get_J (int /*TransitionIndex*/, struct Level */*MyDictionary*/);
__declspec(dllexport) double Partition_Function (double */*Constants*/, double /*Temperature*/);
__declspec(dllexport) double E_tau (int /*TransitionIndex*/, double /*Kappa*/, struct ETauStruct /*ETStruct*/);
__declspec(dllexport) double Rigid_Rotor (double /*A*/, double /*C*/, int /*J*/, int /*Index*/, double /*Kappa*/, struct ETauStruct /*ETStruct*/);
__declspec(dllexport) double Get_Frequency (int /*J_Up*/, int /*J_Low*/, int /*IndexUp*/, int /*IndexLow*/, double */*Constants*/, struct ETauStruct /*ETStruct*/);
__declspec(dllexport) int Get_Catalog (struct Transition */*CatalogtoFill*/, double */*Constants*/, int /*CatLines*/, int /*Verbose*/, struct ETauStruct /*ETStruct*/, struct Level */*MyDictionary*/);

//General catalog functions
__declspec(dllexport) void print_Transition (struct Transition /*TransitionToPrint*/, struct Level */*MyDictionary*/);
__declspec(dllexport) int Fill_Catalog_Restricted (struct Transition */*SourceCatalog*/, struct Transition **/*CatalogtoFill*/, double */*Constants*/, int /*CatLines*/, double /*FrequencyLow*/, double /*FrequencyHigh*/, int */*Dipoles*/, int /*Verbose*/, struct Level */*MyDictionary*/);
__declspec(dllexport) void Sort_Catalog (struct Transition */*Catalog*/, int /*TransitionCount*/, int /*SortMethod*/);
__declspec(dllexport) int Catalog_Comparator (const void */*a*/, const void */*b*/);
__declspec(dllexport) int Catalog_Comparator_Index_Upper (const void */*a*/, const void */*b*/);
__declspec(dllexport) int Catalog_Comparator_Index_Lower (const void */*a*/, const void */*b*/);
__declspec(dllexport) void insertionSort(struct Transition */*CatalogtoSort*/, int /*TransitionCount*/);
__declspec(dllexport) void Calculate_State_Energies (double **/*UpperStateEnergies*/, struct Transition */*SourceCatalog*/, int /*CatalogTransitions*/);
__declspec(dllexport) void Calculate_Intensities (double **/*Intensity*/, struct Transition */*SourceCatalog*/, int /*CatalogTransitions*/, double */*Energies*/, double /*T*/, double */*Dipoles*/);
