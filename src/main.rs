

use eyre::Result;
use rsfibreseq::App;

fn main() -> Result<()> {
    let app = App::new();
    app.run()
}
