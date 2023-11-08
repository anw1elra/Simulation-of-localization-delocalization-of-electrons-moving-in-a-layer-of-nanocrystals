extern crate plotly;
use plotly::color::{NamedColor, Rgb};
use plotly::common::{Anchor, Font, Line, Marker, MarkerSymbol, Mode, Title};
use plotly::histogram::Bins;
use plotly::layout::{Axis, Legend, Shape, ShapeLine, ShapeType, Margin};
use plotly::{ImageFormat, Layout, Plot, Scatter, Histogram};
use std::fs;
use std::path::PathBuf;

fn line_and_scatter_plot(x: Vec<f64>, y: Vec<f64>, flnm: PathBuf, ylab: &str, legend: &str, col: NamedColor) {
    let bgcol = Rgb::new(255, 255, 255);
    let linecol1 = col;
    let forecol = Rgb::new(0, 0, 0);
    let gridcol = Rgb::new(180, 180, 180);
    let transp = NamedColor::Transparent;
    let thick: usize = 3;
    let medium: usize = 3;
    let _thin: usize = 2;
    let msize: usize = 10;
    let fsz_title: usize = 35;
    let fsz_legend: usize = 35;
    let fsz_ticks: usize = 25;
    let fsz_axes: usize = 35;

    let trace1 = Scatter::new(x.clone(), y)
        .name(legend)
        .mode(Mode::LinesMarkers)
        .line(Line::new().color(linecol1).width(medium as f64))
        .marker(Marker::new().size(msize).symbol(MarkerSymbol::Circle));

    let title = Title::new("Electrons vs Quantum dots")
        .font(Font::new().size(fsz_title).family("Serif").color(forecol));

    let legend = Legend::new()
        .x(0.01)
        .x_anchor(Anchor::Left)
        .y(1.0 - 0.0133)
        .y_anchor(Anchor::Top)
        .font(Font::new().size(fsz_legend).color(forecol).family("Serif"))
        .border_width(medium)
        .border_color(forecol)
        .background_color(bgcol)
        .item_width(52);

    let axis = Axis::new()
        .position(0.0)
        .show_line(true)
        .line_color(forecol)
        .line_width(thick)
        .tick_length(6)
        .tick_width(medium)
        .tick_color(forecol)
        .tick_font(Font::new().color(forecol))
        .zero_line(false)
        .show_grid(true)
        .grid_color(gridcol);

    let axisx = axis.clone().title(
        Title::new("Electric field, kV/cm")
            .font(Font::new().size(fsz_axes).color(forecol).family("Serif")));

    let axisy = axis
        .clone()
        .title(Title::new(ylab)
            .font(Font::new().size(fsz_axes).color(forecol).family("Serif")))
            .tick_angle(270.0);

    let line_top = Shape::new()
        .shape_type(ShapeType::Line)
        .x_ref("paper")
        .y_ref("paper")
        .x0(0.)
        .y0(1.)
        .x1(1.)
        .y1(1.)
        .line(ShapeLine::new().color(forecol).width(thick as f64));

    let line_right = Shape::new()
        .shape_type(ShapeType::Line)
        .x_ref("paper")
        .y_ref("paper")
        .x0(1.)
        .y0(0.)
        .x1(1.)
        .y1(1.)
        .line(ShapeLine::new().color(forecol).width(thick as f64));

    let mut layout = Layout::new()
        .width(1024)
        .height(768)
        .font(Font::new().size(fsz_ticks))
        .title(title)
        .legend(legend)
        .show_legend(true)
        .x_axis(axisx)
        .y_axis(axisy)
        .plot_background_color(transp)
        .paper_background_color(bgcol)
        .margin(Margin::new().left(100).bottom(100));

    layout.add_shape(line_top);
    layout.add_shape(line_right);

    let mut plot = Plot::new();
    plot.add_trace(trace1);
    plot.set_layout(layout);

    //plot.write_html(flnm);
    //plot.write_image(flnm, ImageFormat::SVG, 1024, 768, 1.0);
    plot.write_image(flnm, ImageFormat::PNG, 1280, 960, 1.0);
}

fn histogram_plot(y: Vec<f64>, flnm: PathBuf, xlab: &str, title_text: &str, fillcol: NamedColor) {
    let bgcol = Rgb::new(255, 255, 255);
    let forecol = Rgb::new(0, 0, 0);
    let gridcol = Rgb::new(180, 180, 180);
    let transp = NamedColor::Transparent;
    let thick: usize = 3;
    let medium: usize = 3;
    let _thin: usize = 2;
    let fsz_title: usize = 35;
    let fsz_ticks: usize = 25;
    let fsz_axes: usize = 35;

    let bars = Histogram::new(y)
        .name("").x_bins(Bins::new(-0.5, crate::VMAX, 0.025*crate::VMAX))
        .marker(Marker::new().color(fillcol));

    let title = Title::new(title_text)
        .font(Font::new().size(fsz_title).family("Serif").color(forecol));

    let axis = Axis::new()
        .position(0.0)
        .show_line(true)
        .line_color(forecol)
        .line_width(thick)
        .tick_length(6)
        .tick_width(medium)
        .tick_color(forecol)
        .tick_font(Font::new().color(forecol))
        .zero_line(false)
        .show_grid(true)
        .grid_color(gridcol);

    let axisx = axis.clone().title(
        Title::new(xlab)
            .font(Font::new().size(fsz_axes).color(forecol).family("Serif")))
            .fixed_range(true)
            .range(vec![-0.5,crate::VMAX]);

    let axisy = axis
        .clone()
        .title(Title::new("Number of electrons")
            .font(Font::new().size(fsz_axes).color(forecol).family("Serif")))
            .tick_angle(270.0);

    let line_top = Shape::new()
        .shape_type(ShapeType::Line)
        .x_ref("paper")
        .y_ref("paper")
        .x0(0.)
        .y0(1.)
        .x1(1.)
        .y1(1.)
        .line(ShapeLine::new().color(forecol).width(thick as f64));

    let line_right = Shape::new()
        .shape_type(ShapeType::Line)
        .x_ref("paper")
        .y_ref("paper")
        .x0(1.)
        .y0(0.)
        .x1(1.)
        .y1(1.)
        .line(ShapeLine::new().color(forecol).width(thick as f64));

    let mut layout = Layout::new()
        .width(1024)
        .height(768)
        .font(Font::new().size(fsz_ticks))
        .title(title)
        .x_axis(axisx)
        .y_axis(axisy)
        .plot_background_color(transp)
        .paper_background_color(bgcol)
        .margin(Margin::new().left(100).bottom(100))
        .bar_gap(0.1);

    layout.add_shape(line_top);
    layout.add_shape(line_right);

    let mut plot = Plot::new();
    plot.add_trace(bars);
    plot.set_layout(layout);

    //plot.write_html(flnm);
    //plot.write_image(flnm, ImageFormat::SVG, 1024, 768, 1.0);
    plot.write_image(flnm, ImageFormat::PNG, 1280, 960, 1.0);
}

pub fn plot_data() -> std::io::Result<()> {

    //make directory
    let dir_name = "plots";
    fs::create_dir_all(dir_name).expect("Error creating directory");

    let data = crate::file::read_from_file("results","IV_curve.dat");
    let x = data[0].clone();
    let y1 = data[1].clone();
    let y2 = data[2].clone();

    let file_path: PathBuf = [dir_name, "Velocity vs field.png"].iter().collect();
    line_and_scatter_plot(x.clone(), y1, file_path, "Velocity, nm/fs", " <Vx>  ", NamedColor::DarkBlue);
    let file_path: PathBuf = [dir_name, "Energy vs field.png"].iter().collect();
    line_and_scatter_plot(x, y2, file_path, "Energy, eV", " <E>  ", NamedColor::DarkRed);

    for j_field in 1..=crate::IV_POINTS {
        let ee = j_field as f64 * crate::ELECTRIC_MAX / crate::IV_POINTS as f64 * 1e4;
        let fl_name = &format!("Velocity_distribution_E-{}.dat", j_field);
        let data = crate::file::read_from_file("results",fl_name);

        let y1 = data[0].clone();
        let flnm = &format!("Velocity_distribution_E-{}", j_field);
        let file_path: PathBuf = [dir_name, flnm].iter().collect();
        histogram_plot(y1, file_path, "Velocity Vx, nm/fs", &format!("E = {:.3} kV/cm", ee), NamedColor::Chocolate);
    }
    
    Ok(())
}
