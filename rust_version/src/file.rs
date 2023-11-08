pub fn read_from_file(dir_name: &str, fl_name: &str) -> Vec<Vec<f64>> {
    use std::io::BufReader;
    use std::path::PathBuf;
    use std::io::BufRead;
    use std::fs::File;

    let file_path: PathBuf = [dir_name, fl_name].iter().collect();
    let file = File::open(file_path).expect("Error opening file");
    let reader = BufReader::new(file);
    let mut file_contents: Vec<String> = Vec::new();

    for line in reader.lines() {
        file_contents.push(line.unwrap());
    }

    let mut data_rows: Vec<Vec<f64>> = Vec::new();
    let mut data_columns: Vec<Vec<f64>> = Vec::new();

    let mut n_columns: Vec<usize> = Vec::new();
    let mut n_rows: usize = 0;

    for line in file_contents {
        let row: Vec<f64> = line
        .split(|c| c == ' ' || c == '\t'|| c == ',')
        .map(|s| s.trim()) 
        .filter(|s| !s.is_empty()) 
        .map(|s| s.parse().unwrap()) 
        .collect();

        n_columns.push(row.len());
        n_rows = n_rows+1;
        
        data_rows.push(row)
    };

    let n1 = n_columns[1];

    for _i in 0..n1 {
        let column: Vec<f64> = Vec::new();
        data_columns.push(column)
    }

    for j in 0..n_rows {
        for i in 0..n_columns[j] {
            data_columns[i].push(data_rows[j][i])
        }
    }

    data_columns

}